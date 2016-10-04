#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "nb/math_bot.h"
#include "nb/memory_bot.h"
#include "nb/interpolation_bot.h"
#include "nb/eigen_bot.h"
#include "nb/geometric_bot.h"
#include "nb/graph_bot.h"
#include "nb/pde_bot/material.h"
#include "nb/pde_bot/gauss_legendre_quad.h"
#include "nb/pde_bot/common_solid_mechanics/analysis2D.h"
#include "nb/pde_bot/common_solid_mechanics/formulas.h"
#include "nb/pde_bot/boundary_conditions/bcond.h"
#include "nb/pde_bot/boundary_conditions/bcond_iter.h"

#include "set_bconditions.h"

#define INTEGRATOR_TYPE NB_TRIAN

static uint32_t get_cvfa_memsize(uint32_t N_elems);
static void distribute_cvfa_memory(char *memblock,
				   uint32_t N_elems, double **F,
				   nb_partition_t **intmsh);
static void load_integration_mesh(const nb_partition_t *part,
				  nb_partition_t *intmsh);
static void init_global_matrix(const nb_partition_t *intmsh, vcn_sparse_t **K);
static int solver(const vcn_sparse_t *const A,
		  const double *const b, double* x);
static void get_permutation(const vcn_sparse_t *const A,
			    uint32_t *perm, uint32_t *iperm);
static void vector_permutation(uint32_t N, const double *v,
			       const uint32_t *perm, double *vp);

int nb_cvfa_compute_2D_Solid_Mechanics
			(const nb_partition_t *const part,
			 const nb_material_t *const material,
			 const nb_bcond_t *const bcond,
			 bool enable_self_weight, double gravity[2],
			 nb_analysis2D_t analysis2D,
			 nb_analysis2D_params *params2D,
			 double *displacement, /* Output */
			 double *strain,       /* Output */
			 char *boundary_mask   /* Output */)
{
	int status;
	uint32_t N_elems = nb_partition_get_N_elems(part);
	uint32_t memsize = get_cvfa_memsize(N_elems);
	char *memblock = NB_SOFT_MALLOC(memsize);
	double *F;
	nb_partition_t *intmsh;
	distribute_cvfa_memory(memblock, N_elems, &F, &intmsh);

	nb_partition_init(intmsh , INTEGRATOR_TYPE);
	load_integration_mesh(part, intmsh);

	//vcn_sparse_t *K;
	//init_global_matrix(intmsh, &K);

	//nb_graph_init(face_elems_conn);
	//load_face_elems_conn(part, face_elems_conn);
	nb_partition_export_draw(part, "../../../AA_PART.png", 1000, 800,
				 NB_NODE, NB_NULL, NULL, true);/* TEMPORAL */
	nb_partition_export_draw(intmsh, "../../../AA_INTMSH.png", 1000, 800,
				 NB_NODE, NB_NULL, NULL, true);/* TEMPORAL */

	//assemble_global_forces(F, part, material, enable_self_weight,
	//		       gravity);

	//assemble_global_stiffness(K, part, face_elems_conn, material,
	//			  analysis2D, params2D);
	
	//nb_cvfa_set_bconditions(part, material, analysis2D, K, F, bcond, 1.0);

	/*int solver_status = solver(K, F, displacement);
	if (0 != solver_status) {
		status = 1;
		goto CLEANUP_LINEAR_SYSTEM;
	}
	*/

	//compute_strain(strain, boundary_mask, face_elems_conn, part,
	//	       bcond, displacement);

	status = 0;
CLEANUP_LINEAR_SYSTEM:
	//vcn_sparse_destroy(K);
	nb_partition_finish(intmsh);
	NB_SOFT_FREE(memsize, memblock);
	return status;
}

static uint32_t get_cvfa_memsize(uint32_t N_elems)
{
	uint32_t intmsh_size = nb_partition_get_memsize(INTEGRATOR_TYPE);
	uint32_t system_size = 2 * N_elems * sizeof(double);
	return intmsh_size + system_size;
}

static void distribute_cvfa_memory(char *memblock,
				   uint32_t N_elems, double **F,
				   nb_partition_t **intmsh)
{
	uint32_t system_size = 2 * N_elems * sizeof(double);
	*F = (void*) memblock;
	*intmsh = (void*) (memblock + system_size);
}

static void load_integration_mesh(const nb_partition_t *part,
				  nb_partition_t *intmsh)
{
	uint32_t N_elems = nb_partition_get_N_elems(part);

	uint32_t mesh_size = nb_mesh_get_memsize();
	uint32_t vtx_size = 2 * N_elems * sizeof(double);
	uint32_t perm_size = N_elems * sizeof(uint32_t);
	uint32_t memsize = mesh_size + vtx_size + perm_size;
	char *memblock = NB_SOFT_MALLOC(memsize);
	nb_mesh_t *mesh = (void*) memblock;
	double *vtx = (void*) (memblock + mesh_size);
	uint32_t *perm = (void*) (memblock + mesh_size + vtx_size);
	
	for (uint32_t i = 0; i < N_elems; i++) {
		vtx[i * 2] = nb_partition_elem_get_x(part, i);
		vtx[i*2+1] = nb_partition_elem_get_y(part, i);
	}	

	nb_mesh_init(mesh);
	nb_mesh_get_smallest_ns_alpha_complex(mesh, N_elems, vtx, 0.7);
	nb_partition_load_from_mesh(intmsh, mesh);
	nb_mesh_finish(mesh);

	for (uint32_t i = 0; i < N_elems; i++) {
		uint32_t id = nb_partition_get_invtx(intmsh, i);
		perm[id] = i;
	}

	nb_partition_set_nodal_permutation(intmsh, perm);

	NB_SOFT_FREE(memsize, memblock);
}

static void init_global_matrix(const nb_partition_t *intmsh, vcn_sparse_t **K)
{
	uint32_t memsize = nb_graph_get_memsize();
	nb_graph_t *graph = NB_SOFT_MALLOC(memsize);

	nb_graph_init(graph);
	nb_partition_load_graph(intmsh, graph, NB_NODES_LINKED_BY_ELEMS);

	*K = vcn_sparse_create(graph, NULL, 2);

	nb_graph_finish(graph);
	NB_SOFT_FREE(memsize, graph);
}

static int solver(const vcn_sparse_t *const A,
		  const double *const b, double* x)
{
	uint32_t N = vcn_sparse_get_size(A);
	uint32_t memsize = 2 * N * (sizeof(uint32_t) + sizeof(double));
	char *memblock = NB_SOFT_MALLOC(memsize);
	uint32_t *perm = (void*) memblock;
	uint32_t *iperm = (void*) (memblock + N * sizeof(uint32_t));
	double *br = (void*) (memblock + 2 * N * sizeof(uint32_t));
	double *xr = (void*) (memblock + 2 * N * sizeof(uint32_t) +
			      N * sizeof(double));

	double usym = vcn_sparse_get_usym(A);/* TEMPORAL */
	printf("-- K usym: %e\n", usym);     /* TEMPORAL */

	get_permutation(A, perm, iperm);

	vcn_sparse_t *Ar = vcn_sparse_create_permutation(A, perm, iperm);
	vector_permutation(N, b, perm, br);

	int status = vcn_sparse_solve_using_LU(Ar, br, xr, 1);

	vector_permutation(N, xr, iperm, x);
	
	vcn_sparse_destroy(Ar);
	NB_SOFT_FREE(memsize, memblock);
	return status;
}

static void get_permutation(const vcn_sparse_t *const A,
			    uint32_t *perm, uint32_t *iperm)
{
	uint16_t memsize = nb_graph_get_memsize();
	nb_graph_t *graph = NB_SOFT_MALLOC(memsize);
	nb_graph_init(graph);
	nb_sparse_get_graph(A, graph);
	nb_graph_labeling(graph, perm, iperm, NB_LABELING_ND);
	nb_graph_finish(graph);

	NB_SOFT_FREE(memsize, graph);
}

static void vector_permutation(uint32_t N, const double *v,
			       const uint32_t *perm, double *vp)
{
	for (uint32_t i = 0; i < N; i++)
		vp[i] = v[perm[i]];
}

void nb_cvfa_compute_stress_from_strain(const nb_partition_t *part,
					const nb_material_t *const material,
					nb_analysis2D_t analysis2D,
					const double* strain,
					double* stress /* Output */)
{
	uint32_t N_faces = nb_partition_get_N_edges(part);
	for (uint32_t i = 0; i < N_faces; i++) {
		double D[4];
		nb_pde_get_constitutive_matrix(D, material, analysis2D);
		
		stress[i * 3] = (strain[i * 3] * D[0] +
				 strain[i*3+1] * D[1]);
		stress[i*3+1] = (strain[i * 3] * D[1] +
				 strain[i*3+1] * D[2]);
		stress[i*3+2] = strain[i*3+2] * D[3];
	}
}
