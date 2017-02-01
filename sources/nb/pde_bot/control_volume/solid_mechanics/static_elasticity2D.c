#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "nb/math_bot.h"
#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/solver_bot.h"
#include "nb/geometric_bot.h"
#include "nb/graph_bot.h"
#include "nb/pde_bot.h"

#include "../calculation_points.h"
#include "../integration_mesh.h"

#include "elasticity2D.h"
#include "set_bconditions.h"

#define SMOOTH 2

static uint32_t get_cvfa_memsize(uint32_t N_elems, uint32_t N_faces);
static void distribute_cvfa_memory(char *memblock, uint32_t N_elems,
				   uint32_t N_faces, double **xc, double **F,
				   nb_mesh2D_t **intmsh, nb_graph_t **trg_x_vol,
				   face_t ***faces, nb_glquadrature_t *glq);

int nb_cvfa_compute_2D_Solid_Mechanics
			(const nb_mesh2D_t *const mesh,
			 const nb_material_t *const material,
			 const nb_bcond_t *const bcond,
			 bool enable_self_weight, double gravity[2],
			 nb_analysis2D_t analysis2D,
			 nb_analysis2D_params *params2D,
			 double *displacement, /* Output */
			 double *strain,       /* Output */
			 char *boundary_mask   /* Output */)
{
	uint32_t N_elems = nb_mesh2D_get_N_elems(mesh);
	uint32_t N_faces = nb_mesh2D_get_N_edges(mesh);
	uint32_t memsize = get_cvfa_memsize(N_elems, N_faces);
	char *memblock = nb_soft_allocate_mem(memsize);
	double *xc;
	double *F;
	nb_mesh2D_t *intmsh;
	nb_graph_t *trg_x_vol;
	face_t **faces;
	nb_glquadrature_t glq;
	distribute_cvfa_memory(memblock, N_elems, N_faces, &xc, &F,
			       &intmsh, &trg_x_vol, &faces, &glq);

	nb_glquadrature_load(&glq, SMOOTH + 1);

  	nb_cvfa_set_calculation_points(mesh, xc);
	nb_cvfa_init_integration_mesh(intmsh);
	nb_cvfa_load_integration_mesh(intmsh, N_elems, xc, mesh);

	nb_graph_init(trg_x_vol);
	nb_cvfa_correlate_mesh_and_integration_mesh(mesh, intmsh,
						    trg_x_vol);

	nb_sparse_t *K;
	nb_cvfa_init_global_matrix(&K, trg_x_vol, mesh, intmsh, 2);

	nb_cvfa_load_faces(mesh, intmsh, trg_x_vol, faces);

	
	nb_cvfa_assemble_global_forces(F, mesh, material, enable_self_weight,
				       gravity);

	nb_cvfa_assemble_global_stiffness(K, mesh, SMOOTH, intmsh, xc, faces,
					  material, analysis2D, params2D,
					  &glq, NULL);

	nb_cvfa_set_bconditions(mesh, material, analysis2D, 
				K, F, bcond, 1.0);

	int status = nb_sparse_relabel_and_solve_using_LU(K, F,
							  displacement, 1);

	if (status != 0)
		goto CLEAN_AND_EXIT;

	nb_cvfa_compute_strain(strain, boundary_mask, faces, mesh, SMOOTH,
			       intmsh, xc, bcond, displacement, &glq);

	status = 0;
CLEAN_AND_EXIT:
	nb_cvfa_finish_faces(N_faces, faces);
	nb_sparse_destroy(K);
	nb_graph_finish(trg_x_vol);
	nb_mesh2D_finish(intmsh);
	nb_soft_free_mem(memsize, memblock);
	return status;
}

static uint32_t get_cvfa_memsize(uint32_t N_elems, uint32_t N_faces)
{
	uint32_t system_size = 4 * N_elems * sizeof(double);
	uint32_t intmsh_size = nb_cvfa_get_integration_mesh_memsize();
	uint32_t graph_size = nb_graph_get_memsize();
	uint16_t Nq = SMOOTH + 1;
	uint32_t glq_size = 2 * Nq * sizeof(double);
	uint32_t faces_size = N_faces * (sizeof(void*) + sizeof(face_t));
	return graph_size + system_size + intmsh_size + faces_size + glq_size;
}

static void distribute_cvfa_memory(char *memblock, uint32_t N_elems,
				   uint32_t N_faces, double **xc, double **F,
				   nb_mesh2D_t **intmsh, nb_graph_t **trg_x_vol,
				   face_t ***faces, nb_glquadrature_t *glq)
{
	uint32_t system_size = 2 * N_elems * sizeof(double);
	uint32_t intmsh_size = nb_cvfa_get_integration_mesh_memsize();
	uint32_t graph_size = nb_graph_get_memsize();
	uint16_t Nq = SMOOTH + 1;
	uint32_t glq_size = 2 * Nq * sizeof(double);
	*F = (void*) memblock;
	*xc = (void*) (memblock + system_size);
	*intmsh = (void*) (memblock + 2 * system_size);
	*trg_x_vol = (void*) (memblock + 2 * system_size + intmsh_size);
	glq->x = (void*) (memblock + 2 * system_size +
			  intmsh_size + graph_size);
	glq->w = (void*) (memblock + 2 * system_size +
			  intmsh_size + graph_size + Nq * sizeof(double));
	*faces = (void*) (memblock + 2 * system_size +
			  intmsh_size + graph_size + glq_size);
	memblock +=  2 * system_size + intmsh_size + graph_size +
		glq_size + N_faces * sizeof(void*);
	for (uint32_t i = 0; i < N_faces; i++) {
		(*faces)[i] = (void*) (memblock + i * sizeof(face_t));
		memset((*faces)[i], 0, sizeof(face_t));
	}
}

void nb_cvfa_compute_stress_from_strain(const nb_mesh2D_t *mesh,
					const nb_material_t *const material,
					nb_analysis2D_t analysis2D,
					const double* strain,
					double* stress /* Output */)
{
	uint32_t N_faces = nb_mesh2D_get_N_edges(mesh);
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
