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

#include "set_bconditions.h"

#define POW2(a) ((a)*(a))

static int assemble_system(vcn_sparse_t *K, double *F,
			   const nb_partition_t *const part,
			   const nb_material_t *material,
			   bool enable_self_weight,
			   double gravity[2],
			   nb_analysis2D_t analysis2D,
			   nb_analysis2D_params *params2D);
static void assemble_elem(uint32_t elem_id,
			  vcn_sparse_t *K, double *F,
			  const nb_partition_t *const part,
			  const nb_material_t *material,
			  bool enable_self_weight,
			  double gravity[2],
			  nb_analysis2D_t analysis2D,
			  nb_analysis2D_params *params2D);
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
			  vcn_sparse_t *K,  double *F);
static void integrate_inface(uint32_t elem_id, uint16_t face_id,
			     const nb_partition_t *const part,
			     const nb_material_t *material,
			     nb_analysis2D_t analysis2D,
			     nb_analysis2D_params *params2D,
			     vcn_sparse_t *K);
static void get_Ke(const nb_partition_t *const part,
		   uint32_t elem_i, uint16_t ngb_id,
		   const double D[4], double Ke[8]);
static int solver(const vcn_sparse_t *const A,
		  const double *const b, double* x);
static void compute_strain(double *strain,
			   const nb_partition_t *const part,
			   double *disp,
			   nb_analysis2D_t analysis2D,
			   const nb_material_t *const material);
static void estimate_Du(uint32_t elem_id,
			const nb_partition_t *const part,
			double *u, double Du[4]);
static uint16_t get_N_ngb_and_calculate_G(uint32_t elem_id,
					  const nb_partition_t *const part,
					  double *u, double G[4],
					  double Du_sum[4]);

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

	uint32_t N_elems = nb_partition_get_N_elems(part);
	uint32_t F_memsize = 2 * N_elems * sizeof(double);
	double* F = NB_SOFT_MALLOC(F_memsize);
	memset(F, 0, F_memsize);

	int status_assemble = assemble_system(K, F, part, material,
					      enable_self_weight, gravity,
					      analysis2D, params2D);

	if (0 != status_assemble) {
		status = 1;
		goto CLEANUP_LINEAR_SYSTEM;
	}

	nb_cvfa_set_bconditions(part, K, F, bcond, 1.0);

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
			   nb_analysis2D_params *params2D)
{
	uint32_t N_elems = nb_partition_get_N_elems(part);
	vcn_sparse_reset(K);
	memset(F, 0, vcn_sparse_get_size(K) * sizeof(*F));

	for (uint32_t i = 0; i < N_elems; i++) {
		assemble_elem(i, K, F, part, material,
			      enable_self_weight, gravity,
			      analysis2D, params2D);
			
	}

	return 0;
}

static void assemble_elem(uint32_t elem_id,
			  vcn_sparse_t *K, double *F,
			  const nb_partition_t *const part,
			  const nb_material_t *material,
			  bool enable_self_weight,
			  double gravity[2],
			  nb_analysis2D_t analysis2D,
			  nb_analysis2D_params *params2D)
{
	integrate_elem_force(part, material, enable_self_weight,
			     gravity, elem_id, F);

	uint16_t N_faces = nb_partition_elem_get_N_adj(part, elem_id);
	for (uint16_t j = 0; j < N_faces; j++) {
		assemble_face(elem_id, j, part, material,
			      analysis2D, params2D,
			      K, F);
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
			  vcn_sparse_t *K, double *F)
{
	if (nb_partition_elem_has_ngb(part, elem_id, face_id))
		integrate_inface(elem_id, face_id, part, material,
				 analysis2D, params2D, K);
}

static void integrate_inface(uint32_t elem_id, uint16_t face_id,
			     const nb_partition_t *const part,
			     const nb_material_t *material,
			     nb_analysis2D_t analysis2D,
			     nb_analysis2D_params *params2D,
			     vcn_sparse_t *K)
{	
	double D[4];
	nb_pde_get_constitutive_matrix(D, material, analysis2D);

	double Ke[8];
	get_Ke(part, elem_id, face_id, D, Ke);

	double factor = params2D->thickness;

	uint32_t i = elem_id;
	uint32_t j = nb_partition_elem_get_ngb(part, elem_id, face_id);
	vcn_sparse_add(K, i * 2, i * 2, factor * Ke[0]);
	vcn_sparse_add(K, i * 2, i*2+1, factor * Ke[1]);
	vcn_sparse_add(K, i * 2, j * 2, factor * Ke[2]);
	vcn_sparse_add(K, i * 2, j*2+1, factor * Ke[3]);
	vcn_sparse_add(K, i*2+1, i * 2, factor * Ke[4]);
	vcn_sparse_add(K, i*2+1, i*2+1, factor * Ke[5]);
	vcn_sparse_add(K, i*2+1, j * 2, factor * Ke[6]);
	vcn_sparse_add(K, i*2+1, j*2+1, factor * Ke[7]);
}

static void get_Ke(const nb_partition_t *const part,
		   uint32_t elem_i, uint16_t ngb_id,
		   const double D[4], double Ke[8])
{
	double nf[2];
	double lij = nb_partition_elem_face_get_normal(part, elem_i,
						       ngb_id, nf);
	double nij[2];
	double dist = nb_partition_elem_ngb_get_normal(part, elem_i,
						       ngb_id, nij);

	double aij = nij[0] / dist;
	double bij = nij[1] / dist;

	Ke[0] = aij * nf[0] * D[0] + bij * nf[1] * D[3];
	Ke[1] = bij * nf[0] * D[1] + aij * nf[1] * D[3];
	Ke[2] = -Ke[0];
	Ke[3] = -Ke[1];
	Ke[4] = aij * nf[1] * D[1] + bij * nf[0] * D[3];
	Ke[5] = bij * nf[1] * D[2] + aij * nf[0] * D[3];
	Ke[6] = -Ke[4];
	Ke[7] = -Ke[5];

	for (uint8_t k = 0; k < 8; k++)
		Ke[k] *= lij;
}


static int solver(const vcn_sparse_t *const A,
		  const double *const b, double* x)
{
	uint32_t N = vcn_sparse_get_size(A);
	memset(x, 0, N * sizeof(*x));
	uint32_t niter; double tol;
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
	uint32_t N_elems = nb_partition_get_N_elems(part);
	for (uint32_t i = 0; i < N_elems; i++) {
		double Du[4];
		estimate_Du(i, part, disp, Du);
		strain[i * 3] = Du[0];
		strain[i*3+1] = Du[3];
		strain[i*3+2] = (Du[1] + Du[2]);
	}
}

static void estimate_Du(uint32_t elem_id,
			const nb_partition_t *const part,
			double *u, double Du[4])
{
	double G[4], Du_sum[4];
	uint16_t N_ngb = get_N_ngb_and_calculate_G(elem_id, part, u,
						   G, Du_sum);
	
  
	if (1e-5 > fabs(vcn_matrix_2X2_det(G))) {
		/* Only normal information */
		Du[0] = Du_sum[0] / N_ngb;
		Du[1] = Du_sum[1] / N_ngb;
		Du[2] = Du_sum[2] / N_ngb;
		Du[3] = Du_sum[3] / N_ngb;
	} else {
		vcn_matrix_2X2_inverse_destructive(G);      
		Du[0] = G[0] * Du_sum[0] + G[1] * Du_sum[1];
		Du[1] = G[0] * Du_sum[2] + G[1] * Du_sum[3];
		Du[2] = G[2] * Du_sum[0] + G[3] * Du_sum[1];
		Du[3] = G[2] * Du_sum[2] + G[3] * Du_sum[3];
	}
}

static uint16_t get_N_ngb_and_calculate_G(uint32_t elem_id,
					  const nb_partition_t *const part,
					  double *u, double G[4],
					  double Du_sum[4])
{
	uint32_t N_elems = nb_partition_get_N_elems(part);
	memset(G, 0, 4 * sizeof(*G));
	memset(Du_sum, 0, 4 * sizeof(*Du_sum));
	uint16_t N_ngb = 0;
	uint16_t N_adj = nb_partition_elem_get_N_adj(part, elem_id);
	for (uint16_t j = 0; j < N_adj; j++) {
		uint16_t ngb_id = nb_partition_elem_get_ngb(part, elem_id, j);
		if (ngb_id < N_elems) {
			double nij[2];
			double dist =
				nb_partition_elem_ngb_get_normal(part,
								 elem_id,
								 j, nij);

			double uij[2];
			uij[0] = u[ngb_id * 2] - u[elem_id * 2];
			uij[1] = u[ngb_id*2+1] - u[elem_id*2+1];

			G[0] += nij[0] * nij[0];
			G[1] += nij[0] * nij[1];
			G[2] += nij[1] * nij[0];
			G[3] += nij[1] * nij[1];
      
			Du_sum[0] += uij[0] * nij[0] / dist;
			Du_sum[1] += uij[0] * nij[1] / dist;
			Du_sum[2] += uij[1] * nij[0] / dist;
			Du_sum[3] += uij[1] * nij[1] / dist;
      
			N_ngb += 1;
		}
	}
	return N_ngb;
}

void nb_cvfa_compute_stress_from_strain(const nb_partition_t *part,
					const nb_material_t *const material,
					nb_analysis2D_t analysis2D,
					double* strain,
					double* stress /* Output */)
{
	uint32_t N_elems = nb_partition_get_N_elems(part);
	for (uint32_t i = 0; i < N_elems; i++) {
		double D[4];
		nb_pde_get_constitutive_matrix(D, material,
					       analysis2D);

		stress[i * 3] = strain[i * 3] * D[0] +
			strain[i*3+1] * D[1];
		stress[i*3+1] = strain[i * 3] * D[1] +
			strain[i*3+1] * D[2];
		stress[i*3+2] = strain[i*3+2] * D[3];
	}
}
