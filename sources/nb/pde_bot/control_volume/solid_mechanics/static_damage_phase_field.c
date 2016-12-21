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

typedef struct {
	const double *displacement;
	const nb_mesh2D_t *intmsh;
	const int smooth;
} eval_damage_data_t;

static uint32_t get_cvfa_memsize(uint32_t N_elems, uint32_t N_faces);
static void distribute_cvfa_memory(char *memblock, uint32_t N_elems,
				   uint32_t N_faces, double **xc, double **F,
				   double **nodal_damage,
				   nb_mesh2D_t **intmsh, nb_graph_t **trg_x_vol,
				   face_t ***faces, nb_glquadrature_t *glq);
static void init_eval_dmg(nb_cvfa_eval_damage_t * eval_dmg, int smooth,
			  const nb_mesh2D_t *intmsh,
			  const double *displacement);
static double get_damage(const face_t *face, uint16_t subface_id,
			 uint8_t gp, const nb_glquadrature_t *glq,
			 const void *data);
static void get_stress(const nb_material_t *material,
		       nb_analysis2D_t analysis2D,
		       const double strain[3],
		       double stress[3]);
static int solve_damage_equation(const nb_mesh2D_t *mesh,
				 const nb_material_t *material,
				 double *nodal_damage, /* Output */
				 const nb_mesh2D_t *intmsh,
				 const double *xc,
				 face_t **faces,
				 double *H, nb_sparse_t *D,
				 nb_glquadrature_t *glq);
static void compute_damage(double *damage, face_t **faces,
			   const nb_mesh2D_t *const mesh,
			   const nb_mesh2D_t *intmsh, const double *xc,
			   const double *nodal_damage,
			   const double *disp,
			   const nb_glquadrature_t *glq);
static void finish_eval_dmg(nb_cvfa_eval_damage_t * eval_dmg);

int nb_cvfa_compute_2D_damage_phase_field
			(const nb_mesh2D_t *const mesh,
			 const nb_material_t *const material,
			 const nb_bcond_t *const bcond,
			 bool enable_self_weight, double gravity[2],
			 nb_analysis2D_t analysis2D,
			 nb_analysis2D_params *params2D,
			 double *displacement, /* Output */
			 double *strain,       /* Output */
			 double *damage,       /* Output */
			 char *boundary_mask   /* Output */)
{
	uint32_t N_elems = nb_mesh2D_get_N_elems(mesh);
	uint32_t N_faces = nb_mesh2D_get_N_edges(mesh);
	uint32_t memsize = get_cvfa_memsize(N_elems, N_faces);
	char *memblock = nb_soft_allocate_mem(memsize);
	double *xc;
	double *rhs;
	double *nodal_damage;
	nb_mesh2D_t *intmsh;
	nb_graph_t *trg_x_vol;
	face_t **faces;
	nb_glquadrature_t glq;
	distribute_cvfa_memory(memblock, N_elems, N_faces, &xc, &rhs,
			       &nodal_damage, &intmsh, &trg_x_vol,
			       &faces, &glq);

	nb_glquadrature_load(&glq, SMOOTH + 1);

  	nb_cvfa_set_calculation_points(mesh, xc);
	nb_cvfa_init_integration_mesh(intmsh);
	nb_cvfa_load_integration_mesh(intmsh, N_elems, xc);

	nb_graph_init(trg_x_vol);
	nb_cvfa_correlate_mesh_and_integration_mesh(mesh, intmsh,
						    trg_x_vol);
  	nb_sparse_t *A;
	nb_cvfa_init_global_matrix(&A, trg_x_vol, intmsh);

	nb_cvfa_load_faces(mesh, intmsh, trg_x_vol, faces);

	nb_cvfa_eval_damage_t eval_dmg;
	init_eval_dmg(&eval_dmg, SMOOTH, intmsh, displacement);

	nb_cvfa_assemble_global_forces(rhs, mesh, material, enable_self_weight,
				       gravity);

	int status = 0;
	while (1) {
		nb_cvfa_assemble_global_stiffness(A, mesh, SMOOTH, intmsh, xc,
						  faces, material, analysis2D,
						  params2D, &glq, &eval_dmg);
		nb_cvfa_set_bconditions(mesh, material, analysis2D, 
					A, rhs, bcond, 1.0);
		
		status = nb_sparse_relabel_and_solve_using_LU(A, rhs,
							      displacement, 1);

		if (status != 0)
			goto CLEAN_AND_EXIT;

		status = solve_damage_equation(mesh, material, nodal_damage,
					       intmsh, xc, faces,
					       rhs, A, &glq);
		if (status != 0)
			goto CLEAN_AND_EXIT;

		if (/* Residual equal to zero*/1)
			break;
	}
	nb_cvfa_compute_strain(strain, boundary_mask, faces, mesh, SMOOTH,
			       intmsh, xc, bcond, displacement, &glq);
	compute_damage(damage, faces, mesh, intmsh, xc,
		       nodal_damage, displacement, &glq);

CLEAN_AND_EXIT:
	nb_cvfa_finish_faces(N_faces, faces);
	nb_sparse_destroy(A);
	nb_graph_finish(trg_x_vol);
	nb_mesh2D_finish(intmsh);
	nb_soft_free_mem(memsize, memblock);
	finish_eval_dmg(&eval_dmg);
	return status;
}

static uint32_t get_cvfa_memsize(uint32_t N_elems, uint32_t N_faces)
{
	uint32_t system_size = 5 * N_elems * sizeof(double);
	uint32_t intmsh_size = nb_cvfa_get_integration_mesh_memsize();
	uint32_t graph_size = nb_graph_get_memsize();
	uint16_t Nq = SMOOTH + 1;
	uint32_t glq_size = 2 * Nq * sizeof(double);
	uint32_t faces_size = N_faces * (sizeof(void*) + sizeof(face_t));
	return graph_size + system_size + intmsh_size + faces_size + glq_size;
}

static void distribute_cvfa_memory(char *memblock, uint32_t N_elems,
				   uint32_t N_faces, double **xc, double **rhs,
				   double **nodal_damage,
				   nb_mesh2D_t **intmsh, nb_graph_t **trg_x_vol,
				   face_t ***faces, nb_glquadrature_t *glq)
{
	uint32_t elem_size = N_elems * sizeof(double);
	uint32_t system_size = 2 * elem_size;
	uint32_t intmsh_size = nb_cvfa_get_integration_mesh_memsize();
	uint32_t graph_size = nb_graph_get_memsize();
	uint16_t Nq = SMOOTH + 1;
	uint32_t glq_size = 2 * Nq * sizeof(double);
	*rhs = (void*) memblock;
	*xc = (void*) (memblock + system_size);
	*nodal_damage = (void*) (memblock + 2 * system_size);
	*intmsh = (void*) (memblock + 2 * system_size + elem_size);
	*trg_x_vol = (void*) (memblock + 2 * system_size +
			      elem_size + intmsh_size);
	glq->x = (void*) (memblock + 2 * system_size + elem_size +
			  intmsh_size + graph_size);
	glq->w = (void*) (memblock + 2 * system_size + elem_size +
			  intmsh_size + graph_size + Nq * sizeof(double));
	*faces = (void*) (memblock + 2 * system_size + elem_size +
			  intmsh_size + graph_size + glq_size);
	memblock +=  2 * system_size + elem_size + intmsh_size + graph_size +
		glq_size + N_faces * sizeof(void*);
	for (uint32_t i = 0; i < N_faces; i++) {
		(*faces)[i] = (void*) (memblock + i * sizeof(face_t));
		memset((*faces)[i], 0, sizeof(face_t));
	}
}

static void init_eval_dmg(nb_cvfa_eval_damage_t * eval_dmg, int smooth,
			  const nb_mesh2D_t *intmsh,
			  const double *displacement)
{
	eval_damage_data_t *data = nb_allocate_mem(sizeof(*data));
	data->displacement = displacement;
	data->intmsh = intmsh;
	eval_dmg->data = data;
	eval_dmg->get_damage = get_damage;
}

static double get_damage(const face_t *face, uint16_t subface_id,
			 uint8_t gp, const nb_glquadrature_t *glq,
			 const void *data)
{
	
	double strain[3];
	memset(strain, 0, 3 * sizeof(*strain));
	nb_cvfa_subface_sum_strain(smooth, intmsh, face, subface, xc,
				   disp, glq, gq, strain);
	double stress[3];
	get_stress(material, analysis2D, strain, stress);

	double vm_stress = nb_pde_get_vm_stress(stress[0], stress[1], stress[2]);
	return (vm_stress>1)?1:vm_stress;
}

static void get_stress(const nb_material_t *material,
		       nb_analysis2D_t analysis2D,
		       const double strain[3],
		       double stress[3])
{
	double D[4];
	nb_pde_get_constitutive_matrix(D, material, analysis2D);
	
	stress[0] = (strain[0] * D[0] + strain[1] * D[1]);
	stress[1] = (strain[0] * D[1] + strain[1] * D[2]);
	stress[2] =  strain[2] * D[3];
}

static int solve_damage_equation(const nb_mesh2D_t *mesh,
				 const nb_material_t *material,
				 double *nodal_damage, /* Output */
				 const nb_mesh2D_t *intmsh,
				 const double *xc,
				 face_t **faces,
				 double *H, nb_sparse_t *D,
				 nb_glquadrature_t *glq)
{
	return 0;/* TEMPORAL */
}

static void compute_damage(double *damage, face_t **faces,
			   const nb_mesh2D_t *const mesh,
			   const nb_mesh2D_t *intmsh, const double *xc,
			   const double *nodal_damage,
			   const double *disp,
			   const nb_glquadrature_t *glq)
{
	/* TEMPORAL */
}

static void finish_eval_dmg(nb_cvfa_eval_damage_t * eval_dmg)
{
	nb_free_mem(eval_dmg->data);
}
