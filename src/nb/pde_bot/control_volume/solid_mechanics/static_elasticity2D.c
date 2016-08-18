#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "nb/memory_bot.h"
#include "nb/eigen_bot.h"
#include "nb/geometric_bot.h"
#include "nb/graph_bot.h"
#include "nb/pde_bot/material.h"
#include "nb/pde_bot/common_solid_mechanics/analysis2D.h"
#include "nb/pde_bot/common_solid_mechanics/formulas.h"
#include "nb/pde_bot/boundary_conditions/bcond.h"
#include "nb/pde_bot/boundary_conditions/bcond_iter.h"

#define POW2(a) ((a)*(a))

static int assemble_system(vcn_sparse_t *K, double *F,
			   const nb_partition_t *const part,
			   const nb_material_t *material,
			   bool enable_self_weight,
			   double gravity[2],
			   nb_analysis2D_t analysis2D,
			   nb_analysis2D_params *params2D,
			   nb_bcond_t *bcond);
static void set_numeric_bconditions(vcn_sparse_t *K, double *F,
				    const nb_partition_t *const part,
				    nb_bcond_t *bcond);

static void set_numeric_bcond_dirichlet(const nb_partition_t *part,
					vcn_sparse_t* K, uint8_t N_dof,
					double* F,
					nb_bcond_iter_t *iter,
					uint32_t elem_id);
static void assemble_elem(uint32_t elem_id,
			  vcn_sparse_t *K, double *F,
			  const nb_partition_t *const part,
			  const nb_material_t *material,
			  bool enable_self_weight,
			  double gravity[2],
			  nb_analysis2D_t analysis2D,
			  nb_analysis2D_params *params2D,
			  nb_bcond_t *bcond,
			  bool *have_dirichlet);

static void integrate_elem_force(const nb_partition_t *part,
				 const nb_material_t *material,
				 bool enable_self_weight,
				 double gravity[2],
				 uint32_t elem_id,
				 double *F);

static void assemble_face(uint32_t elem_id, uint16_t face_id,
			  const nb_partition_t *const part,
			  const nb_material_t *material,
			  nb_analysis2D_t analysis2D,
			  nb_analysis2D_params *params2D,
			  vcn_sparse_t *K,  double *F,
			  nb_bcond_t *bcond,
			  bool is_dirichlet[2]);

static void integrate_inface(uint32_t elem_id, uint16_t face_id,
			     const nb_partition_t *const part,
			     const nb_material_t *material,
			     nb_analysis2D_t analysis2D,
			     nb_analysis2D_params *params2D,
			     vcn_sparse_t *K);

static void get_Ke(const double D[4], const double xi[2],
		   const double xj[2], double Ke[8]);

static int set_boundary_conditions(const nb_partition_t *const part,
				   vcn_sparse_t *K, double *F,
				   const nb_bcond_t *const bcond,
				   double factor);
static void get_bcond_force_on_face(const nb_bcond_t *bcond,
				    const nb_partition_t *part,
				    uint32_t elem_id, uint16_t face_id,
				    double force[2], bool is_dirichlet[2]);
static void get_force_neumann_on_sgm(const nb_bcond_t *bcond,
				     const nb_partition_t *part,
				     uint32_t elem_id, uint16_t face_id,
				     double force[2]);
static void get_force_neumann_on_vtx(const nb_bcond_t *bcond,
				     const nb_partition_t *part,
				     uint32_t elem_id, uint16_t face_id,
				     double force[2]);
static void get_force_dirichlet_on_sgm(const nb_bcond_t *bcond,
				       const nb_partition_t *part,
				       uint32_t elem_id, uint16_t face_id,
				       double force[2],
				       bool is_dirichlet[2]);
static void get_force_dirichlet_on_vtx(const nb_bcond_t *bcond,
				       const nb_partition_t *part,
				       uint32_t elem_id, uint16_t face_id,
				       double force[2],
				       bool is_dirichlet[2]);
static int solver(const vcn_sparse_t *const A,
		  const double *const b, double* x);

static void compute_strain(double *strain,
			   const nb_partition_t *const part,
			   double *disp,
			   nb_analysis2D_t analysis2D,
			   const nb_material_t *const material);

int nb_cvfa_compute_2D_Solid_Mechanics
			(const nb_partition_t *const part,
			 const nb_material_t *const material,
			 const nb_bcond_t *const bcond,
			 bool enable_self_weight,
			 double gravity[2],
			 nb_analysis2D_t analysis2D,
			 nb_analysis2D_params *params2D,
			 double *displacement, /* Output */
			 double *strain       /* Output */)
{
	int status = 0;
	vcn_graph_t *graph = malloc(nb_graph_get_memsize());
	nb_graph_init(graph);
	nb_partition_load_interelem_graph(part, graph);
	vcn_sparse_t *K = vcn_sparse_create(graph, NULL, 2);
	nb_graph_finish(graph);

	uint32_t N_nod = nb_partition_get_N_nodes(part);
	uint32_t F_memsize = 2 * N_nod * sizeof(double);
	double* F = NB_SOFT_MALLOC(F_memsize);
	memset(F, 0, F_memsize);

	int status_assemble = assemble_system(K, F, part, material,
					      enable_self_weight, gravity,
					      analysis2D, params2D, bcond);
	if (0 != status_assemble) {
		status = 1;
		goto CLEANUP_LINEAR_SYSTEM;
	}
  
	int solver_status = solver(K, F, displacement);
	if (0 != solver_status) {
		status = 2;
		goto CLEANUP_LINEAR_SYSTEM;
	}

	compute_strain(strain, part, displacement,
		       analysis2D, material);

	status = 0;
CLEANUP_LINEAR_SYSTEM:
	vcn_sparse_destroy(K);
	NB_SOFT_FREE(F_memsize, F);
	return status;
}

static int assemble_system(vcn_sparse_t *K, double *F,
			   const nb_partition_t *const part,
			   const nb_material_t *material,
			   bool enable_self_weight,
			   double gravity[2],
			   nb_analysis2D_t analysis2D,
			   nb_analysis2D_params *params2D,
			   nb_bcond_t *bcond)
{
	uint32_t N_elems = nb_partition_get_N_elems(part);
	vcn_sparse_reset(K);
	memset(F, 0, vcn_sparse_get_size(K) * sizeof(*F));

	uint16_t bcond_size = nb_bcond_get_memsize(2);
	nb_bcond_t *numeric_bcond = NB_SOFT_MALLOC(bcond_size);
	nb_bcond_init(bcond, 2);
	for (uint32_t i = 0; i < N_elems; i++) {
		bool have_dirichlet[2];
		assemble_elem(i, K, F, part, material,
			      enable_self_weight, gravity,
			      analysis2D, params2D, bcond,
			      have_dirichlet);
		if (have_dirichlet[0] || have_dirichlet[1]) {
			nb_bcond_push(numeric_bcond,
				      NB_DIRICHLET,
				      NB_BC_ON_POINT,
				      i, have_dirichlet,
				      &(F[i * 2]));
		}
			
	}

	set_numeric_bconditions(K, F, part, numeric_bcond);

	NB_SOFT_FREE(bcond_size, numeric_bcond);
	return 0;
}

static void set_numeric_bconditions(vcn_sparse_t *K, double *F,
				    const nb_partition_t *const part,
				    nb_bcond_t *bcond)
{
	uint8_t N_dof = nb_bcond_get_N_dof(bcond);
	uint16_t size = nb_bcond_iter_get_memsize();
	nb_bcond_iter_t *iter = alloca(size);
	nb_bcond_iter_init(iter);
	nb_bcond_iter_set_conditions(iter, bcond,
				     NB_DIRICHLET,
				     NB_BC_ON_POINT);
	while (nb_bcond_iter_has_more(iter)) {
		nb_bcond_iter_go_next(iter);
		uint32_t elem_id = nb_bcond_iter_get_id(iter);
		set_numeric_bcond_dirichlet(part, K, N_dof,
					    F, iter, elem_id);
	}
	nb_bcond_iter_finish(iter);
}

static void set_numeric_bcond_dirichlet(const nb_partition_t *part,
					vcn_sparse_t* K, uint8_t N_dof,
					double* F,
					nb_bcond_iter_t *iter,
					uint32_t elem_id)
{
	double x[2];
	x[0] = nb_partition_get_x_elem(part, elem_id);
	x[1] = nb_partition_get_y_elem(part, elem_id);

	double *val = alloca(N_dof * sizeof(*val));
	nb_bcond_iter_get_val(iter, N_dof, x, 0, val);
	for (uint8_t j = 0; j < N_dof; j++) {
		bool mask = nb_bcond_iter_get_mask(iter, j);
		if (mask) {
			uint32_t mtx_id = elem_id * N_dof + j;
			vcn_sparse_set_Dirichlet_condition(K, F, mtx_id,
							   val[j]);
		}
	}
}

static void assemble_elem(uint32_t elem_id,
			  vcn_sparse_t *K, double *F,
			  const nb_partition_t *const part,
			  const nb_material_t *material,
			  bool enable_self_weight,
			  double gravity[2],
			  nb_analysis2D_t analysis2D,
			  nb_analysis2D_params *params2D,
			  nb_bcond_t *bcond,
			  bool *have_dirichlet)
{
	integrate_elem_force(part, material, enable_self_weight,
			     gravity, i, F);

	have_dirichlet[0] = false;
	have_dirichlet[1] = false;
	uint16_t N_faces = nb_partition_elem_get_N_adj(part, i);
	for (uint16_t j = 0; j < N_faces; j++) {
		bool is_dirichlet[2];
		assemble_face(i, j, part, material,
			      analysis2D, params2D,
			      K, F, bcond, is_dirichlet);

		if (is_dirichlet[0])
			have_dirichlet[0] = true;
		if (is_dirichlet[1])
			have_dirichlet[1] = true;
	}
}

static void integrate_elem_force(const nb_partition_t *part,
				 const nb_material_t *material,
				 bool enable_self_weight,
				 double gravity[2],
				 uint32_t elem_id,
				 double *F)
{
	if (enable_self_weight) {
		double area = nb_partition_elem_get_area(part, elem_id);
		double mass = area * nb_material_get_density(material);
		F[elem_id * 2] += mass * gravity[0];
		F[elem_id*2+1] += mass * gravity[1];
	}
}

static void assemble_face(uint32_t elem_id, uint16_t face_id,
			  const nb_partition_t *const part,
			  const nb_material_t *material,
			  nb_analysis2D_t analysis2D,
			  nb_analysis2D_params *params2D,
			  vcn_sparse_t *K, double *F,
			  nb_bcond_t *bcond,
			  bool is_dirichlet[2])
{
	if (nb_partition_elem_has_ngb(part, elem_id, face_id)) {
		integrate_inface(elem_id, face_id, part, material,
				 analysis2D, params2D, K);
	} else {
		double fij[2];
		get_bcond_force_on_face(bcond, part, elem_id, face_id,
					fij, is_dirichlet);
		F[elem_id * 2] -= fij[0];
		F[elem_id*2+1] -= fij[1];
	}
}

static void integrate_inface(uint32_t elem_id, uint16_t face_id,
			     const nb_partition_t *const part,
			     const nb_material_t *material,
			     nb_analysis2D_t analysis2D,
			     nb_analysis2D_params *params2D,
			     vcn_sparse_t *K)
{
	uint32_t i = elem_id;
	uint32_t j = nb_partition_elem_get_ngb(part, elem_id, face_id);

	double xi[2];
	xi[0] = nb_partition_get_x_elem(part, i);
	xi[1] = nb_partition_get_y_elem(part, i);
	double xj[2];
	xj[0] = nb_partition_get_x_elem(part, j);
	xj[1] = nb_partition_get_y_elem(part, j);
	
	double D[4];
	nb_pde_get_constitutive_matrix(D, material, analysis2D);

	double Ke[8];
	get_Ke(D, xi, xj, Ke);

	double lij = nb_partition_elem_face_get_length(part, elem_id,
						       face_id);
	double factor = lij * params2D->thickness;

	vcn_sparse_add(K, i * 2, i * 2, factor * Ke[0]);
	vcn_sparse_add(K, i * 2, i*2+1, factor * Ke[1]);
	vcn_sparse_add(K, i * 2, j * 2, factor * Ke[2]);
	vcn_sparse_add(K, i * 2, j*2+1, factor * Ke[3]);
	vcn_sparse_add(K, i*2+1, i * 2, factor * Ke[4]);
	vcn_sparse_add(K, i*2+1, i*2+1, factor * Ke[5]);
	vcn_sparse_add(K, i*2+1, j * 2, factor * Ke[6]);
	vcn_sparse_add(K, i*2+1, j*2+1, factor * Ke[7]);
}

static void get_Ke(const double D[4], const double xi[2],
		   const double xj[2], double Ke[8])
{
	double dist2 = POW2(xi[0] - xj[0]) + POW2(xi[1] - xj[1]);
	double aij = (xi[0] - xj[0]) / dist2;
	double bij = (xi[1] - xj[1]) / dist2;

	double dist = sqrt(dist2);
	double n[2];
	n[0] = (xj[0] - xi[0]) / dist;
	n[1] = (xj[1] - xi[1]) / dist;

	Ke[0] = aij * n[0] * D[0] + bij * n[1] * D[3];
	Ke[1] = bij * n[0] * D[1] + aij * n[1] * D[3];
	Ke[2] = -Ke[0];
	Ke[3] = -Ke[1];
	Ke[4] = aij * n[1] * D[1] + bij * n[0] * D[3];
	Ke[5] = bij * n[1] * D[2] + aij * n[0] * D[3];
	Ke[6] = -Ke[4];
	Ke[7] = -Ke[5];
}

static void get_bcond_force_on_face(const nb_bcond_t *bcond,
				    const nb_partition_t *part,
				    uint32_t elem_id, uint16_t face_id,
				    double force[2],
				    bool is_dirichlet[2])
{
	memset(force, 0, 2 * sizeof(*force));
	get_force_neumann_on_sgm(bcond, part, elem_id, face_id, force);
	get_force_neumann_on_vtx(bcond, part, elem_id, face_id, force);
	get_force_dirichlet_on_sgm(bcond, part, elem_id, face_id, force,
				   is_dirichlet);
	get_force_dirichlet_on_vtx(bcond, part, elem_id, face_id, force,
				   is_dirichlet);
}

static void get_force_neumann_on_sgm(const nb_bcond_t *bcond,
				     const nb_partition_t *part,
				     uint32_t elem_id, uint16_t face_id,
				     double force[2])
{
	;/* PENDING */
}

static void get_force_neumann_on_vtx(const nb_bcond_t *bcond,
				     const nb_partition_t *part,
				     uint32_t elem_id, uint16_t face_id,
				     double force[2])
{
	;/* PENDING */
}

static void get_force_dirichlet_on_sgm(const nb_bcond_t *bcond,
				       const nb_partition_t *part,
				       uint32_t elem_id, uint16_t face_id,
				       double force[2],
				       bool is_dirichlet[2])
{
	;/* PENDING */
}

static void get_force_dirichlet_on_vtx(const nb_bcond_t *bcond,
				       const nb_partition_t *part,
				       uint32_t elem_id, uint16_t face_id,
				       double force[2],
				       bool is_dirichlet[2])
{
	;/* PENDING y NPI de como hacerlo */
}

static int solver(const vcn_sparse_t *const A,
		  const double *const b, double* x)
{
	uint32_t N = vcn_sparse_get_size(A);
	memset(x, 0, N * sizeof(*x));
	int status = vcn_sparse_solve_CG_precond_Jacobi(A, b, x, N,
							1e-8, NULL,
							NULL, 1);
	int out;
	if (0 == status || 1 == status)
		out = 0;
	else
		out = 1; /* Tolerance not reached in CG Jacobi */
	return out;
}

static void compute_strain(double *strain,
			   const nb_partition_t *const part,
			   double *disp,
			   nb_analysis2D_t analysis2D,
			   const nb_material_t *const material)
{
	;/* PENDING */
}
