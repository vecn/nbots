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
#include "nb/pde_bot/common_solid_mechanics/boundary_conditions.h"
#include "nb/pde_bot/boundary_conditions/bcond.h"
#include "nb/pde_bot/boundary_conditions/bcond_iter.h"

#define POW2(a) ((a)*(a))

static int assemble_system(vcn_sparse_t *K, double *F,
			   const nb_mshpoly_t *const mesh,
			   const nb_material_t *material,
			   bool enable_self_weight,
			   double gravity[2],
			   nb_analysis2D_t analysis2D,
			   nb_analysis2D_params *params2D);

static void assemble_face(uint32_t elem_id, uint16_t face_id,
			  const nb_mshpoly_t *const mesh,
			  const nb_material_t *material,
			  bool enable_self_weight,
			  double gravity[2],
			  nb_analysis2D_t analysis2D,
			  nb_analysis2D_params *params2D,
			  vcn_sparse_t *K, double *F);
static void get_Ke(const double D[4], const double xi[2],
		   const double xj[2], double Ke[8]);
static double get_face_length(uint32_t elem_id, uint16_t face_id,
			      const nb_mshpoly_t *const mesh);
static int set_boundary_conditions(const nb_mshpoly_t *const mesh,
				   vcn_sparse_t *K, double *F,
				   const nb_bcond_t *const bcond,
				   double factor);
static int solver(const vcn_sparse_t *const A,
		  const double *const b, double* x);

static void compute_strain(double *strain,
			   const nb_mshpoly_t *const mesh,
			   double *disp,
			   nb_analysis2D_t analysis2D,
			   const nb_material_t *const material);

int nb_cvfa_compute_2D_Solid_Mechanics
			(const nb_mshpoly_t *const mesh,
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
	nb_mshpoly_set_elemental_graph(mesh, graph);
	vcn_sparse_t *K = vcn_sparse_create(graph, NULL, 2);
	nb_graph_finish(graph);

	uint32_t F_memsize = 2 * mesh->N_nod * sizeof(double);
	double* F = NB_SOFT_MALLOC(F_memsize);
	memset(F, 0, F_memsize);

	int status_assemble = assemble_system(K, F, mesh, material,
					      enable_self_weight, gravity,
					      analysis2D, params2D);
	if (0 != status_assemble) {
		status = 1;
		goto CLEANUP_LINEAR_SYSTEM;
	}

	set_boundary_conditions(mesh, K, F, bcond, 1.0);

  
	int solver_status = solver(K, F, displacement);
	if (0 != solver_status) {
		status = 2;
		goto CLEANUP_LINEAR_SYSTEM;
	}

	compute_strain(strain, mesh, displacement,
		       analysis2D, material);

	status = 0;
CLEANUP_LINEAR_SYSTEM:
	vcn_sparse_destroy(K);
	NB_SOFT_FREE(F_memsize, F);
	return status;
}

static int assemble_system(vcn_sparse_t *K, double *F,
			   const nb_mshpoly_t *const mesh,
			   const nb_material_t *material,
			   bool enable_self_weight,
			   double gravity[2],
			   nb_analysis2D_t analysis2D,
			   nb_analysis2D_params *params2D)
{
	uint32_t N_elems = mesh->N_elems;
	vcn_sparse_reset(K);
	memset(F, 0, vcn_sparse_get_size(K) * sizeof(*F));

	for (uint32_t i = 0; i < N_elems; i++) {
		uint16_t N_faces = mesh->N_adj[i];
		for (uint16_t j = 0; j < N_faces; j++) {
			assemble_face(i, j, mesh, material,
				      enable_self_weight, gravity,
				      analysis2D, params2D, K, F);
		}
	}
	return 0;
}

static void assemble_face(uint32_t elem_id, uint16_t face_id,
			  const nb_mshpoly_t *const mesh,
			  const nb_material_t *material,
			  bool enable_self_weight,
			  double gravity[2],
			  nb_analysis2D_t analysis2D,
			  nb_analysis2D_params *params2D,
			  vcn_sparse_t *K, double *F)
{
	uint32_t i = elem_id;
	uint32_t j = mesh->ngb[elem_id][face_id];
	double *xi = &(mesh->cen[i * 2]);
	double *xj = &(mesh->cen[j * 2]);

	double D[4];
	nb_pde_get_constitutive_matrix(D, material, analysis2D);

	double Ke[8];
	get_Ke(D, xi, xj, Ke);

	double lij = get_face_length(elem_id, face_id, mesh);

	vcn_sparse_add(K, i * 2, i * 2, Ke[0]);
	vcn_sparse_add(K, i * 2, i*2+1, Ke[1]);
	vcn_sparse_add(K, i * 2, j * 2, Ke[2]);
	vcn_sparse_add(K, i * 2, j*2+1, Ke[3]);
	vcn_sparse_add(K, i*2+1, i * 2, Ke[4]);
	vcn_sparse_add(K, i*2+1, i*2+1, Ke[5]);
	vcn_sparse_add(K, i*2+1, j * 2, Ke[6]);
	vcn_sparse_add(K, i*2+1, j*2+1, Ke[7]);
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

static double get_face_length(uint32_t elem_id, uint16_t face_id,
			      const nb_mshpoly_t *const mesh)
{
	uint32_t id1 = mesh->adj[elem_id][face_id];
	uint16_t next_node = (face_id + 1) % mesh->N_adj[elem_id];
	uint32_t id2 = mesh->adj[elem_id][next_node];

	double *v1 = &(mesh->nod[id1*2]);
	double *v2 = &(mesh->nod[id2*2]);

	return sqrt(POW2(v2[0] - v1[0]) + POW2(v2[1] - v1[1]));
}

static int set_boundary_conditions(const nb_mshpoly_t *const mesh,
				   vcn_sparse_t *K, double *F,
				   const nb_bcond_t *const bcond,
				   double factor)
{
	return 1;
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
			   const nb_mshpoly_t *const mesh,
			   double *disp,
			   nb_analysis2D_t analysis2D,
			   const nb_material_t *const material)
{
	;
}
