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

#define MAX_ITER 100
#define SMOOTH 1

#define POW2(a) ((a)*(a))

enum {
	ELASTIC_SOLVER_FAILS = 1,
	DAMAGE_SOLVER_FAILS,
	NO_INVERSE
};

typedef struct {
	int smooth;
	const nb_mesh2D_t *intmsh;
	const double *disp;
	const double *xc;
	const nb_material_t *material;
	nb_analysis2D_t analysis2D;
} eval_damage_data_t;

static uint32_t get_cvfa_memsize(uint32_t N_elems, uint32_t N_faces);
static void distribute_cvfa_memory(char *memblock, uint32_t N_elems,
				   uint32_t N_faces, double **xc, double **F,
				   double **delta_disp, double **residual,
				   double **nodal_damage, double **rhs_damage,
				   nb_mesh2D_t **intmsh, nb_graph_t **trg_x_vol,
				   face_t ***faces, nb_glquadrature_t *glq);
static void init_eval_dmg(nb_cvfa_eval_damage_t * eval_dmg, int smooth,
			  const nb_mesh2D_t *intmsh,
			  const double *displacement,
			  const double *xc,
			  const nb_material_t *material,
			  nb_analysis2D_t analysis2D);
static double get_damage(const face_t *face, uint16_t subface_id,
			 uint8_t gp, const nb_glquadrature_t *glq,
			 const void *data);
static int minimize_residual(const nb_mesh2D_t *const mesh,
			     const nb_material_t *const material,
			     const nb_bcond_t *const bcond,
			     bool enable_self_weight, double gravity[2],
			     nb_analysis2D_t analysis2D,
			     nb_analysis2D_params *params2D,
			     double *displacement, /* Output */
			     double *strain,       /* Output */
			     double *elem_damage,  /* Output */
			     const nb_mesh2D_t *intmsh,
			     const double *xc,
			     face_t **faces,
			     nb_sparse_t *A,
			     nb_glquadrature_t *glq,
			     const nb_cvfa_eval_damage_t *eval_dmg,
			     double *F,
			     double *delta_disp,
			     double *residual,
			     double *rhs_damage,
			     double bc_factor,
			     uint32_t id_elem_monitor,
			     double *monitor_reaction);
static int solve_damage_equation(const nb_mesh2D_t *mesh,
				 const nb_material_t *material,
				 double *nodal_damage, /* Output */
				 const nb_mesh2D_t *intmsh,
				 const double *xc,
				 face_t **faces,
				 double *H, nb_sparse_t *D,
				 nb_glquadrature_t *glq);

static void show_error_message(int status);
static void save_reaction_log(const char *logfile, uint32_t iter,
			      double factor, double reaction);
static void compute_damage(double *damage, face_t **faces,
			   const nb_mesh2D_t *const mesh,
			   const double *elem_damage,
			   const nb_glquadrature_t *glq,
			   const nb_cvfa_eval_damage_t *eval_dmg);
static void get_face_damage(double *damage,
			    face_t **faces, uint32_t face_id,
			    const nb_mesh2D_t *mesh,
			    const double *elem_damage,
			    const nb_glquadrature_t *glq,
			    const nb_cvfa_eval_damage_t *eval_dmg);
static void get_internal_face_damage(double *damage,
				     face_t **faces, uint32_t face_id,
				     const nb_glquadrature_t *glq,
				     const nb_cvfa_eval_damage_t *eval_dmg);
static void get_boundary_face_damage(double *damage, face_t **faces,
				     uint32_t face_id);
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
	double *F;
	double *delta_disp;
	double *residual;
	double *rhs_damage;
	double *elem_damage;
	nb_mesh2D_t *intmsh;
	nb_graph_t *trg_x_vol;
	face_t **faces;
	nb_glquadrature_t glq;
	distribute_cvfa_memory(memblock, N_elems, N_faces, &xc, &F,
			       &delta_disp, &residual, &elem_damage, 
			       &rhs_damage, &intmsh, &trg_x_vol,
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
	init_eval_dmg(&eval_dmg, SMOOTH, intmsh, displacement, xc,
		      material, analysis2D);

	uint32_t id_elem_monitor[2];
	nb_cvfa_get_elem_adj_to_model_node(mesh, 8, id_elem_monitor);

	int status = 0;
	memset(displacement, 0, 2 * N_elems * sizeof(*displacement));

	uint32_t N_steps = 5;
	double step_factor = 1.0 / N_steps;
	for (int i = 0; i < N_steps; i++) {
		double bc_factor = (i+1) * step_factor;
		double reaction;
		status = minimize_residual(mesh, material, bcond,
					   enable_self_weight,
					   gravity, analysis2D, params2D,
					   displacement, strain, elem_damage,
					   intmsh, xc, faces, A, &glq,
					   &eval_dmg, F, delta_disp, residual,
					   rhs_damage, bc_factor,
					   id_elem_monitor[0], &reaction);
		if (status != 0) {
			show_error_message(status);
			goto CLEAN_AND_EXIT;
		}

		save_reaction_log("Reaction.log", i, bc_factor, reaction);
	}

	nb_cvfa_compute_strain(strain, boundary_mask, faces, mesh, SMOOTH,
			       intmsh, xc, bcond, displacement, &glq);
	compute_damage(damage, faces, mesh, elem_damage, &glq, &eval_dmg);

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
	uint32_t system_size = 10 * N_elems * sizeof(double);
	uint32_t intmsh_size = nb_cvfa_get_integration_mesh_memsize();
	uint32_t graph_size = nb_graph_get_memsize();
	uint16_t Nq = SMOOTH + 1;
	uint32_t glq_size = 2 * Nq * sizeof(double);
	uint32_t faces_size = N_faces * (sizeof(void*) + sizeof(face_t));
	return graph_size + system_size + intmsh_size + faces_size + glq_size;
}

static void distribute_cvfa_memory(char *memblock, uint32_t N_elems,
				   uint32_t N_faces, double **xc, double **F,
				   double **delta_disp, double **residual,
				   double **elem_damage, double **rhs_damage,
				   nb_mesh2D_t **intmsh, nb_graph_t **trg_x_vol,
				   face_t ***faces, nb_glquadrature_t *glq)
{
	uint32_t elem_size = N_elems * sizeof(double);
	uint32_t system_size = 2 * elem_size;
	uint32_t intmsh_size = nb_cvfa_get_integration_mesh_memsize();
	uint32_t graph_size = nb_graph_get_memsize();
	uint16_t Nq = SMOOTH + 1;
	uint32_t glq_size = 2 * Nq * sizeof(double);
	*F = (void*) memblock;
	*xc = (void*) (memblock + system_size);
	*delta_disp = (void*) (memblock + 2 * system_size);
	*residual = (void*) (memblock + 3 * system_size);
	*elem_damage = (void*) (memblock + 4 * system_size);
	*rhs_damage = (void*) (memblock + 4 * system_size + elem_size);
	*intmsh = (void*) (memblock + 4 * system_size + 2 * elem_size);
	*trg_x_vol = (void*) (memblock + 4 * system_size +
			      2 * elem_size + intmsh_size);
	glq->x = (void*) (memblock + 4 * system_size + 2 * elem_size +
			  intmsh_size + graph_size);
	glq->w = (void*) (memblock + 4 * system_size + 2 * elem_size +
			  intmsh_size + graph_size + Nq * sizeof(double));
	*faces = (void*) (memblock + 4 * system_size + 2 * elem_size +
			  intmsh_size + graph_size + glq_size);
	memblock +=  4 * system_size + 2 * elem_size + intmsh_size +
		graph_size + glq_size + N_faces * sizeof(void*);
	for (uint32_t i = 0; i < N_faces; i++) {
		(*faces)[i] = (void*) (memblock + i * sizeof(face_t));
		memset((*faces)[i], 0, sizeof(face_t));
	}
}

static void init_eval_dmg(nb_cvfa_eval_damage_t * eval_dmg, int smooth,
			  const nb_mesh2D_t *intmsh,
			  const double *displacement,
			  const double *xc,
			  const nb_material_t *material,
			  nb_analysis2D_t analysis2D)
{
	eval_damage_data_t *data = nb_allocate_mem(sizeof(*data));
	data->smooth = smooth;
	data->intmsh = intmsh;
	data->disp = displacement;
	data->xc = xc;
	data->material = material;
	data->analysis2D = analysis2D;
	eval_dmg->data = data;
	eval_dmg->get_damage = get_damage;
}

#define MIN(a,b) (((a)<(b))?(a):(b)) /* TEMPORAL */
static double get_damage(const face_t *face, uint16_t subface_id,
			 uint8_t gp, const nb_glquadrature_t *glq,
			 const void *data)
{
	const eval_damage_data_t* dmg_data = data;

	double strain[3];
	nb_cvfa_subface_get_strain(dmg_data->smooth,
				   dmg_data->intmsh,
				   face, face->subfaces[subface_id],
				   dmg_data->xc,
				   dmg_data->disp,
				   glq, gp, strain);

	double lame[2];
	nb_pde_get_lame_params(lame, dmg_data->material,
			       dmg_data->analysis2D);
	double tr = strain[0] + strain[1];
	double norm2 = POW2(strain[0]) +
		2 * POW2(0.5 * strain[2]) + POW2(strain[1]);
	double energy = lame[0] * norm2 + 0.5 * lame[1] * POW2(tr);

	double h = nb_material_get_damage_length_scale(dmg_data->material);
	double G = nb_material_get_fracture_energy(dmg_data->material);
	return MIN(1.0, h * energy / G);
}

static int minimize_residual(const nb_mesh2D_t *const mesh,
			     const nb_material_t *const material,
			     const nb_bcond_t *const bcond,
			     bool enable_self_weight, double gravity[2],
			     nb_analysis2D_t analysis2D,
			     nb_analysis2D_params *params2D,
			     double *displacement, /* Output */
			     double *strain,       /* Output */
			     double *elem_damage,  /* Output */
			     const nb_mesh2D_t *intmsh,
			     const double *xc,
			     face_t **faces,
			     nb_sparse_t *A,
			     nb_glquadrature_t *glq,
			     const nb_cvfa_eval_damage_t *eval_dmg,
			     double *F,
			     double *delta_disp,
			     double *residual,
			     double *rhs_damage,
			     double bc_factor,
			     uint32_t id_elem_monitor,
			     double *monitor_reaction)
{
	uint16_t bcond_size = nb_bcond_get_memsize(2);
	nb_bcond_t *numeric_bcond = nb_soft_allocate_mem(bcond_size);
	nb_bcond_init(numeric_bcond, 2);
	nb_cvfa_get_numeric_bconditions(mesh, material, analysis2D, bcond,
					bc_factor, numeric_bcond);

	uint32_t N_elems = nb_mesh2D_get_N_elems(mesh);
	int status = 0;
	int iter = 0;
	double rnorm = 1;
	while (iter < MAX_ITER) {
		nb_cvfa_assemble_global_forces(F, mesh, material,
					       enable_self_weight,
					       gravity);
		nb_cvfa_assemble_global_stiffness(A, mesh, SMOOTH, intmsh, xc,
						  faces, material, analysis2D,
						  params2D, glq, eval_dmg);

		nb_cvfa_set_numeric_bconditions(A, F, mesh, numeric_bcond);
		
		nb_sparse_multiply_vector(A, displacement, residual, 1);
		nb_vector_substract_to(2 * N_elems, residual, F);
		rnorm = nb_vector_get_norm(residual, 2 * N_elems);
		if (rnorm < 1e-6)
			goto GET_REACTION;

		if (rnorm != rnorm) {
			status = NO_INVERSE;
			goto EXIT;
		}

		status = nb_sparse_relabel_and_solve_using_LU(A, residual,
							      delta_disp, 2);
		if (status != 0) {
			status = ELASTIC_SOLVER_FAILS;
			goto EXIT;
		}

		nb_vector_sum(2 * N_elems, displacement, delta_disp);

		status = solve_damage_equation(mesh, material, elem_damage,
					       intmsh, xc, faces,
					       rhs_damage, A, glq);
		if (status != 0) {
			status = DAMAGE_SOLVER_FAILS;
			goto EXIT;
		}
		
		iter ++;
		printf(" ----> DAMAGE ITER: %i (%e)\n", iter, rnorm);/* TEMP */
	}
GET_REACTION:
	nb_sparse_multiply_vector(A, displacement, residual, 1);
	*monitor_reaction = nb_math_hypo(residual[id_elem_monitor * 2],
					 residual[id_elem_monitor*2+1]);
EXIT:
	printf(" >>>>> [%i] DAMAGE ITER: %i (%e)\n",
	       status, iter, rnorm);/* TEMPORAL */

	nb_bcond_finish(numeric_bcond);
	nb_soft_free_mem(bcond_size, numeric_bcond);
	return status;
}

static int solve_damage_equation(const nb_mesh2D_t *mesh,
				 const nb_material_t *material,
				 double *elem_damage, /* Output */
				 const nb_mesh2D_t *intmsh,
				 const double *xc,
				 face_t **faces,
				 double *H, nb_sparse_t *D,
				 nb_glquadrature_t *glq)
{
	return 0;/* TEMPORAL */
}

static void show_error_message(int status)
{
	fprintf(stderr, " NB >> Sorry, something goes wrong :(\n");
	fprintf(stderr, " NB >> at PDE_BOT::CONTROL_VOLUME::"	\
		"SOLID_MECHANICS::DAMAGE_PHASE_FIELD\n");
	fprintf(stderr, " NB >> Description: ");
	switch(status) {
	case ELASTIC_SOLVER_FAILS:
		fprintf(stderr, "Elastic solver fails");
		break;
	case DAMAGE_SOLVER_FAILS:
		fprintf(stderr, "Damage solver fails");
		break;
	case NO_INVERSE:
		fprintf(stderr, "Stiffness matrix is singular");
		break;
	default:
		fprintf(stderr, "Unkown error");		
	}
	fprintf(stderr, "\n");
}

static void save_reaction_log(const char *logfile, uint32_t iter,
			      double factor, double reaction)
{
	FILE *fp;
	if (0 == iter)
		fp = fopen(logfile, "w");
	else
		fp = fopen(logfile, "a");

	if (NULL == fp)
		goto EXIT;

	fprintf(fp, "%e %e\n", factor, reaction);
	fclose(fp);
EXIT:
	return;
}

static void compute_damage(double *damage, face_t **faces,
			   const nb_mesh2D_t *const mesh,
			   const double *elem_damage,
			   const nb_glquadrature_t *glq,
			   const nb_cvfa_eval_damage_t *eval_dmg)
{
	uint32_t N_faces = nb_mesh2D_get_N_edges(mesh);

 	for (uint32_t i = 0; i < N_faces; i++)
		get_face_damage(damage, faces, i, mesh,
				elem_damage, glq, eval_dmg);
}

static void get_face_damage(double *damage,
			    face_t **faces, uint32_t face_id,
			    const nb_mesh2D_t *mesh,
			    const double *elem_damage,
			    const nb_glquadrature_t *glq,
			    const nb_cvfa_eval_damage_t *eval_dmg)
{
	uint32_t N_elems = nb_mesh2D_get_N_elems(mesh);
	if (faces[face_id]->elems[1] < N_elems)
		get_internal_face_damage(damage, faces, face_id,
					 glq, eval_dmg);
	else
		get_boundary_face_damage(damage, faces, face_id);
}

static void get_internal_face_damage(double *damage,
				     face_t **faces, uint32_t face_id,
				     const nb_glquadrature_t *glq,
				     const nb_cvfa_eval_damage_t *eval_dmg)
{
	damage[face_id] = 0.0;
	face_t *face = faces[face_id];
	for (uint16_t i = 0; i < face->N_sf; i++) {
		subface_t *subface = face->subfaces[i];
		double lf = nb_utils2D_get_dist(subface->x1, subface->x2);
		for (uint8_t q = 0; q < glq->N; q++) {
			double dmg = eval_dmg->get_damage(face, i, q, glq,
							  eval_dmg->data);
			double wq = lf * glq->w[q] * 0.5;
			damage[face_id] += wq * dmg;
		}
	}
	double length = nb_utils2D_get_dist(face->x1, face->x2);
	damage[face_id] /= length;
}

static void get_boundary_face_damage(double *damage, face_t **faces,
				     uint32_t face_id)
{
	damage[face_id] = 0.0;
}

static void finish_eval_dmg(nb_cvfa_eval_damage_t * eval_dmg)
{
	nb_free_mem(eval_dmg->data);
}
