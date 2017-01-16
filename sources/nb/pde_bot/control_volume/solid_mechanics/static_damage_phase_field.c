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

#define SMOOTH 0
#define MIDPOINT_VOL_INTEGRALS false

#define POW2(a) ((a)*(a))
#define MIN(a,b) (((a)<(b))?(a):(b))

enum {
	SUCCESS_RESIDUAL_MIN = 0,
	ELASTIC_SOLVER_FAILS,
	DAMAGE_SOLVER_FAILS,
	NO_INVERSE,
	MAX_ITER_REACHED,
	LU_ALLOCATION_FAILS,
	LU_ALLOCATION_FAILS_DMG
};

typedef struct {
	int smooth;
	const nb_mesh2D_t *mesh;
	const nb_mesh2D_t *intmsh;
	const double *disp;
	const double *xc;
	const double *elem_dmg;
	const nb_material_t *material;
	nb_analysis2D_t analysis2D;
} eval_damage_data_t;

static uint32_t get_cvfa_memsize(uint32_t N_elems, uint32_t N_faces);
static void distribute_cvfa_memory(char *memblock, uint32_t N_elems,
				   uint32_t N_faces, double **xc,
				   double **nodal_damage,
				   nb_mesh2D_t **intmsh, nb_graph_t **trg_x_vol,
				   face_t ***faces, nb_glquadrature_t *glq,
				   char **minimize_residuals_memblock);
static void init_eval_dmg(nb_cvfa_eval_damage_t * eval_dmg, int smooth,
			  const nb_mesh2D_t *mesh,
			  const nb_mesh2D_t *intmsh,
			  const double *displacement,
			  const double *xc,
			  const double *elem_dmg,
			  const nb_material_t *material,
			  nb_analysis2D_t analysis2D);
static double get_damage(const face_t *face, uint16_t subface_id,
			 uint8_t gp, const nb_glquadrature_t *glq,
			 const void *data);
static double get_internal_subface_damage(const face_t *face,
					  uint16_t subface_id,
					  uint8_t gp,
					  const nb_glquadrature_t *glq,
					  const eval_damage_data_t *dmg_data);
static double get_subface_energy(const face_t *face, uint16_t subface_id,
				 const double xq[2],
				 const eval_damage_data_t *dmg_data);
static double get_elem_energy(uint32_t elem_id,
			      const eval_damage_data_t *dmg_data);
static double get_energy(const double strain[3],
			 const nb_material_t *material,
			 nb_analysis2D_t analysis2D);
static double subface_get_damage_simplexwise
				(const face_t *face,
				 uint16_t subface_id,
				 uint8_t gp,
				 const nb_glquadrature_t *glq,
				 const eval_damage_data_t *dmg_data);
static double subface_get_damage_pairwise(const face_t *face,
					  uint16_t subface_id,
					  uint8_t gp,
					  const nb_glquadrature_t *glq,
					  const eval_damage_data_t *dmg_data);
static uint32_t get_memsize_for_minimize_residual(uint32_t N_elems);
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
			     face_t **faces, int smooth,
			     nb_sparse_t *K, nb_sparse_t *D,
			     const nb_glquadrature_t *glq,
			     const nb_cvfa_eval_damage_t *eval_dmg,
			     double bc_factor, uint32_t max_iter,
			     uint32_t id_elem_monitor,
			     double *monitor_reaction,
			     char *memblock);
static int solve_linear_system(const nb_sparse_t *A,
			       const uint32_t *perm,
			       const uint32_t *iperm,
			       nb_sparse_t *Ar,
			       nb_sparse_t *Lr,
			       nb_sparse_t *Ur,
			       double *b, double *x);
static void assemble_global_damage(const nb_cvfa_eval_damage_t * eval_dmg,
				   const nb_glquadrature_t *glq,
				   face_t **faces, nb_sparse_t *D, double *H);
static void assemble_elem_damage(const eval_damage_data_t *dmg_data,
				 uint32_t elem_id,
				 const nb_glquadrature_t *glq,
				 nb_sparse_t *D, double *H);
static void assemble_face_damage(const eval_damage_data_t *dmg_data,
				 face_t *face,
				 const nb_glquadrature_t *glq,
				 nb_sparse_t *D, double *H);
static void integrate_subvolume_damage(const eval_damage_data_t *dmg_data,
				       face_t *face, uint16_t subface_id,
				       const nb_glquadrature_t *glq,
				       nb_sparse_t *D, double *H);
static void integrate_lateral_subvolume_damage
				(const double st1[2],
				 const double st2[2],
				 const double st3[2],
				 const eval_damage_data_t *dmg_data,
				 const face_t *face, uint32_t subface_id,
				 uint32_t elem_id,
				 const nb_glquadrature_t *glq,
				 nb_sparse_t *D, double *H);

static void integrate_svol_simplexwise_gp_dmg
				(const eval_damage_data_t *dmg_data,
				 const face_t *face, uint32_t subface_id,
				 uint32_t elem_id, double wq,
				 const double xq[2],
				 nb_sparse_t *D, double *H);
static void integrate_svol_pairwise_gp_dmg
				(const eval_damage_data_t *dmg_data,
				 const face_t *face, uint32_t subface_id,
				 uint32_t elem_id, double wq,
				 const double xq[2],
				 nb_sparse_t *D, double *H);
static void integrate_subface_damage(const eval_damage_data_t *dmg_data,
				     face_t *face, uint16_t subface_id,
				     const nb_glquadrature_t *glq,
				     nb_sparse_t *D, double *H);
static void integrate_subface_simplexwise_gp_damage
					(const eval_damage_data_t *dmg_data,
					 face_t *face, uint16_t subface_id,
					 const nb_glquadrature_t *glq,
					 uint8_t gp,
					 nb_sparse_t *D, double *H);
static void integrate_subface_pairwise_gp_damage
					(const eval_damage_data_t *dmg_data,
					 face_t *face,
					 uint16_t subface_id,
					 const nb_glquadrature_t *glq,
					 uint8_t gp,
					 nb_sparse_t *D, double *H);
static void integrate_exterior_subface_damage
					(const eval_damage_data_t *dmg_data,
					 face_t *face,
					 const nb_glquadrature_t *glq,
					 nb_sparse_t *D, double *H);
static void show_error_message(int status);
static void get_reaction_log_name(const char *dir, char *reaction_log);
static void save_reaction_log(const char *logfile, uint32_t iter,
			      double factor, double reaction);
static void save_simulation(const char *dir,
			    const nb_mesh2D_t *mesh, const double *disp,
			    const double *elem_damage,
			    const double *face_damage, int iter, double time);
static void save_data(const char *name,
		      const nb_mesh2D_t *mesh,
		      const double *disp,
		      const double *elem_damage,
		      const double *face_damage, int iter, double time,
		      const char *mesh_name);
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
			 const char *dir_to_save,
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
	double *elem_damage;
	nb_mesh2D_t *intmsh;
	nb_graph_t *trg_x_vol;
	face_t **faces;
	nb_glquadrature_t glq;
	char *minimize_residual_memblock;
	distribute_cvfa_memory(memblock, N_elems, N_faces, &xc,
			       &elem_damage, &intmsh, &trg_x_vol,
			       &faces, &glq, &minimize_residual_memblock);

	nb_glquadrature_load(&glq, SMOOTH + 1);

  	nb_cvfa_set_calculation_points(mesh, xc);
	nb_cvfa_init_integration_mesh(intmsh);
	nb_cvfa_load_integration_mesh(intmsh, N_elems, xc);

	nb_graph_init(trg_x_vol);
	nb_cvfa_correlate_mesh_and_integration_mesh(mesh, intmsh,
						    trg_x_vol);
  	nb_sparse_t *K;
	nb_cvfa_init_global_matrix(&K, trg_x_vol, intmsh, 2);

	nb_cvfa_load_faces(mesh, intmsh, trg_x_vol, faces);

	nb_sparse_t *D;
	nb_cvfa_init_global_matrix(&D, trg_x_vol, intmsh, 1);

	nb_cvfa_eval_damage_t eval_dmg;
	init_eval_dmg(&eval_dmg, SMOOTH, mesh, intmsh, displacement,
		      xc, elem_damage, material, analysis2D);

	uint32_t id_elem_monitor[2];
	nb_cvfa_get_elem_adj_to_model_node(mesh, 8, id_elem_monitor);

	int status = 0;
	memset(displacement, 0, 2 * N_elems * sizeof(*displacement));
	memset(elem_damage, 0, N_elems * sizeof(*elem_damage));

	char reaction_log[100];
	get_reaction_log_name(dir_to_save, reaction_log);
	save_reaction_log(reaction_log, 0, 0, 0);
	uint32_t iter = 1;
	double bc_factor_increment = 0.25;
	uint32_t max_iter = 20;
	double bc_factor = 0;
	while (bc_factor < 1.0) {
		bc_factor += bc_factor_increment;
		double reaction;
		double trim_bc_factor = MIN(1.0, bc_factor);
		status = minimize_residual(mesh, material, bcond,
					   enable_self_weight,
					   gravity, analysis2D, params2D,
					   displacement, strain, elem_damage,
					   intmsh, xc, faces, SMOOTH, K, D,
					   &glq, &eval_dmg, trim_bc_factor,
					   max_iter,
					   id_elem_monitor[0], &reaction,
					   minimize_residual_memblock);
		if (status == MAX_ITER_REACHED) {
			bc_factor -= bc_factor_increment;
			bc_factor_increment /= 2;
			max_iter *= 2;
		} else {
			if (status != 0) {
				show_error_message(status);
				goto CLEAN_AND_EXIT;
			}

			save_reaction_log(reaction_log, iter,
					  MIN(1.0, bc_factor),
					  reaction);
			if (max_iter > 25) {
				bc_factor_increment *= 2;
				max_iter /= 2;
			}
			compute_damage(damage, faces, mesh, elem_damage,
				       &glq, &eval_dmg);
			save_simulation(dir_to_save, mesh, displacement,
					elem_damage, damage, iter,
					trim_bc_factor);
			iter ++;
		}
	}

	nb_cvfa_compute_strain(strain, boundary_mask, faces, mesh, SMOOTH,
			       intmsh, xc, bcond, displacement, &glq);

CLEAN_AND_EXIT:
	nb_cvfa_finish_faces(N_faces, faces);
	nb_sparse_destroy(K);
	nb_sparse_destroy(D);
	nb_graph_finish(trg_x_vol);
	nb_mesh2D_finish(intmsh);
	nb_soft_free_mem(memsize, memblock);
	finish_eval_dmg(&eval_dmg);
	return status;
}

static uint32_t get_cvfa_memsize(uint32_t N_elems, uint32_t N_faces)
{
	uint32_t system_size = 3 * N_elems * sizeof(double);
	uint32_t intmsh_size = nb_cvfa_get_integration_mesh_memsize();
	uint32_t graph_size = nb_graph_get_memsize();
	uint16_t Nq = SMOOTH + 1;
	uint32_t glq_size = 2 * Nq * sizeof(double);
	uint32_t faces_size = N_faces * (sizeof(void*) + sizeof(face_t));
	uint32_t minimize_residual = get_memsize_for_minimize_residual(N_elems);
	return graph_size + system_size + intmsh_size + faces_size +
		glq_size + minimize_residual;
}

static void distribute_cvfa_memory(char *memblock, uint32_t N_elems,
				   uint32_t N_faces, double **xc,
				   double **elem_damage, nb_mesh2D_t **intmsh,
				   nb_graph_t **trg_x_vol,
				   face_t ***faces, nb_glquadrature_t *glq,
				   char **minimize_residuals_memblock)
{
	uint32_t elem_size = N_elems * sizeof(double);
	uint32_t system_size = 2 * elem_size;
	uint32_t intmsh_size = nb_cvfa_get_integration_mesh_memsize();
	uint32_t graph_size = nb_graph_get_memsize();
	uint16_t Nq = SMOOTH + 1;
	uint32_t glq_size = 2 * Nq * sizeof(double);
	*xc = (void*) memblock;
	*elem_damage = (void*) (memblock + system_size);
	*intmsh = (void*) (memblock + system_size + elem_size);
	*trg_x_vol = (void*) (memblock + system_size +
			      elem_size + intmsh_size);
	glq->x = (void*) (memblock + system_size + elem_size +
			  intmsh_size + graph_size);
	glq->w = (void*) (memblock + system_size + elem_size +
			  intmsh_size + graph_size + Nq * sizeof(double));
	*faces = (void*) (memblock + system_size + elem_size +
			  intmsh_size + graph_size + glq_size);
	memblock +=  system_size + elem_size + intmsh_size +
		graph_size + glq_size + N_faces * sizeof(void*);
	for (uint32_t i = 0; i < N_faces; i++) {
		(*faces)[i] = (void*) memblock;
		memblock += sizeof(face_t);
		memset((*faces)[i], 0, sizeof(face_t));
	}
	*minimize_residuals_memblock = memblock;
}

static void init_eval_dmg(nb_cvfa_eval_damage_t * eval_dmg, int smooth,
			  const nb_mesh2D_t *mesh,
			  const nb_mesh2D_t *intmsh,
			  const double *displacement,
			  const double *xc,
			  const double *elem_dmg,
			  const nb_material_t *material,
			  nb_analysis2D_t analysis2D)
{
	eval_damage_data_t *data = nb_allocate_mem(sizeof(*data));
	data->smooth = smooth;
	data->mesh = mesh;
	data->intmsh = intmsh;
	data->disp = displacement;
	data->xc = xc;
	data->elem_dmg = elem_dmg;
	data->material = material;
	data->analysis2D = analysis2D;
	eval_dmg->data = data;
	eval_dmg->get_damage = get_damage;
}

static double get_damage(const face_t *face, uint16_t subface_id,
			 uint8_t gp, const nb_glquadrature_t *glq,
			 const void *data)
{
	const eval_damage_data_t *dmg_data = data;

	double damage = 0.0;
	if (nb_cvfa_face_is_internal(face, dmg_data->mesh)) {
		damage = get_internal_subface_damage(face, subface_id,
						     gp, glq, dmg_data);
	}
	return damage;
}

static double get_internal_subface_damage(const face_t *face,
					  uint16_t subface_id,
					  uint8_t gp,
					  const nb_glquadrature_t *glq,
					  const eval_damage_data_t *dmg_data)
{
	const subface_t *subface = face->subfaces[subface_id];
	double damage;
	if (nb_cvfa_subface_in_simplex(subface))
		damage = subface_get_damage_simplexwise(face, subface_id, gp,
							glq, dmg_data);
	else
		damage = subface_get_damage_pairwise(face, subface_id, gp,
						     glq, dmg_data);
	return damage;
}

static double subface_get_damage_simplexwise(const face_t *face,
					     uint16_t subface_id,
					     uint8_t gp,
					     const nb_glquadrature_t *glq,
					     const eval_damage_data_t *dmg_data)
{
	const subface_t *subface = face->subfaces[subface_id];

	double t1[2], t2[2], t3[2];
	nb_cvfa_load_trg_points(dmg_data->intmsh,
				subface->trg_id, t1, t2, t3);

	double xq[2];
	nb_cvfa_subface_get_xq(subface, glq, gp, xq);
	double xi[2];
	nb_cvfa_get_normalized_point(dmg_data->smooth, t1, t2, t3, xq, xi);

	double damage = 0.0;
	for (uint8_t k = 0; k < 3; k++) {
		uint32_t elem_id =
			nb_mesh2D_elem_get_adj(dmg_data->intmsh, 
					       subface->trg_id, k);
		double xk;
		if (k < 2)
			xk = xi[k];
		else
			xk = 1.0 - xi[0] - xi[1];
		damage += xk * dmg_data->elem_dmg[elem_id];
	}
	return damage;
}

static double get_subface_energy(const face_t *face, uint16_t subface_id,
				 const double xq[2],
				 const eval_damage_data_t *dmg_data)
{
	double strain[3];
	nb_cvfa_subface_get_strain(dmg_data->smooth,
				   dmg_data->intmsh,
				   face, face->subfaces[subface_id],
				   dmg_data->xc,
				   dmg_data->disp,
				   xq, strain);
	
	return get_energy(strain, dmg_data->material, dmg_data->analysis2D);
}

static double get_elem_energy(uint32_t elem_id,
			      const eval_damage_data_t *dmg_data)
{
	const nb_mesh2D_t *mesh = dmg_data->mesh;
	double x[2];
	x[0] = dmg_data->xc[elem_id * 2];
	x[1] = dmg_data->xc[elem_id*2+1];
	double u[2];
	u[0] = dmg_data->disp[elem_id * 2];
	u[1] = dmg_data->disp[elem_id*2+1];

	uint16_t N_adj = nb_mesh2D_elem_get_N_adj(mesh, elem_id);
	double ni[30];/* WARNING: Max 15 neighbours */
	double ui[30];
	uint16_t N = 0;
	for (uint32_t i = 0; i < N_adj; i++) {
		if (nb_mesh2D_elem_has_ngb(mesh, elem_id, i)) {
			uint32_t ngb =
				nb_mesh2D_elem_get_ngb(mesh, elem_id, i);
			ni[N * 2] = dmg_data->xc[ngb * 2];
			ni[N*2+1] = dmg_data->xc[ngb*2+1];
			ui[N * 2] = dmg_data->disp[ngb * 2];
			ui[N*2+1] = dmg_data->disp[ngb*2+1];
			N ++;
		}		
	}

	double Du[4];
	nb_pde_get_frechet_derivative(N, 2, 2, x, u, ni, ui, Du);

	double strain[3];
	strain[0] = Du[0];
	strain[1] = Du[3];
	strain[2] = Du[1] + Du[2];
	
	return get_energy(strain, dmg_data->material, dmg_data->analysis2D);	
}

static double get_energy(const double strain[3],
			 const nb_material_t *material,
			 nb_analysis2D_t analysis2D)
{
	double lame[2];
	nb_pde_get_lame_params(lame, material, analysis2D);
	double tr = strain[0] + strain[1];
	double norm2 = POW2(strain[0]) +
		2 * POW2(0.5 * strain[2]) + POW2(strain[1]);
	double energy = lame[0] * norm2 + 0.5 * lame[1] * POW2(tr);
	return energy;
}

static double subface_get_damage_pairwise(const face_t *face,
					  uint16_t subface_id,
					  uint8_t gp,
					  const nb_glquadrature_t *glq,
					  const eval_damage_data_t *dmg_data)
{
	const subface_t *subface = face->subfaces[subface_id];

	uint32_t id1 = face->elems[0];
	uint32_t id2 = face->elems[1];
	double c1[2], c2[2];
	c1[0] = dmg_data->xc[id1 * 2];
	c1[1] = dmg_data->xc[id1*2+1];
	c2[0] = dmg_data->xc[id2 * 2];
	c2[1] = dmg_data->xc[id2*2+1];

	double xq[2];
	nb_cvfa_subface_get_xq(subface, glq, gp, xq);

	double xdiff = c2[0] - c1[0];
	double ydiff = c2[1] - c1[1];
	double d2 = nb_utils2D_get_dist2(c1, c2);
	double dot = (xq[0] - c1[0]) * xdiff + (xq[1] - c1[1]) * ydiff;
	double z = dot / d2;

	double damage = 0.0;
	for (uint8_t k = 0; k < 2; k++) {
		uint32_t elem_id = face->elems[k];
		if (k == 1)
			z = 1 - z;
		damage += z * dmg_data->elem_dmg[elem_id];
	}
	return damage;
}

static uint32_t get_memsize_for_minimize_residual(uint32_t N_elems)
{
	uint32_t memsize = 6 * N_elems * sizeof(uint32_t) +
		9 * N_elems * sizeof(double);
	return memsize;
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
			     face_t **faces, int smooth,
			     nb_sparse_t *K, nb_sparse_t *D,
			     const nb_glquadrature_t *glq,
			     const nb_cvfa_eval_damage_t *eval_dmg,
			     double bc_factor, uint32_t max_iter,
			     uint32_t id_elem_monitor,
			     double *monitor_reaction,
			     char *memblock)
{
	int status = SUCCESS_RESIDUAL_MIN;
	uint32_t N = nb_sparse_get_size(K);
	uint32_t *perm = (void*) memblock;
	uint32_t *iperm = (void*) (memblock + N * sizeof(uint32_t));
	double *F = (void*) (memblock + 2 * N * sizeof(uint32_t));
	double *delta_disp = (void*) (memblock + 2 * N * sizeof(uint32_t) +
				      N * sizeof(double));
	double *residual = (void*) (memblock + 2 * N * sizeof(uint32_t) +
				    2 * N * sizeof(double));
	double *aux_disp = (void*) (memblock + 2 * N * sizeof(uint32_t) +
				    3 * N * sizeof(double));
	uint32_t N_elems = nb_mesh2D_get_N_elems(mesh);
	uint32_t *dmg_perm = (void*) (memblock + 2 * N * sizeof(uint32_t) +
				      4 * N * sizeof(double));
	uint32_t *dmg_iperm = (void*) (memblock + 2 * N * sizeof(uint32_t) +
				       4 * N * sizeof(double) +
				       N_elems * sizeof(uint32_t));
	double *rhs_damage = (void*) (memblock + 2 * N * sizeof(uint32_t) +
				      4 * N * sizeof(double) +
				      2 * N_elems * sizeof(uint32_t));

	uint16_t bcond_size = nb_bcond_get_memsize(2);
	nb_bcond_t *numeric_bcond = nb_soft_allocate_mem(bcond_size);
	nb_bcond_init(numeric_bcond, 2);
	nb_cvfa_get_numeric_bconditions(mesh, material, analysis2D, bcond,
					bc_factor, numeric_bcond);

	memcpy(aux_disp, displacement, N * sizeof(*aux_disp));

	nb_sparse_calculate_permutation(K, perm, iperm);

	nb_sparse_t *Kr = nb_sparse_create_permutation(K, perm, iperm);
	nb_sparse_t *Lr = NULL; 
	nb_sparse_t *Ur = NULL;
	nb_sparse_alloc_LU(Kr, &Lr, &Ur);
	if(NULL == Lr) {
		status = LU_ALLOCATION_FAILS;
		goto EXIT;
	}

	nb_sparse_calculate_permutation(D, dmg_perm, dmg_iperm);

	nb_sparse_t *Dr = nb_sparse_create_permutation(D, dmg_perm, dmg_iperm);
	nb_sparse_t *dmg_Lr = NULL; 
	nb_sparse_t *dmg_Ur = NULL;
	nb_sparse_alloc_LU(Dr, &dmg_Lr, &dmg_Ur);
	if(NULL == dmg_Lr) {
		status = LU_ALLOCATION_FAILS_DMG;
		goto EXIT;
	}

	int iter = 0;
	double rnorm = 1;
	while (1) {
		nb_cvfa_assemble_global_forces(F, mesh, material,
					       enable_self_weight,
					       gravity);
		nb_cvfa_assemble_global_stiffness(K, mesh, smooth,
						  intmsh, xc, faces, material,
						  analysis2D, params2D, glq,
						  eval_dmg);

		nb_cvfa_set_numeric_bconditions(K, F, mesh, numeric_bcond);
		
		nb_sparse_multiply_vector(K, displacement, residual, 1);
		nb_vector_substract_to(2 * N_elems, residual, F);
		rnorm = nb_vector_get_norm(residual, 2 * N_elems);
		if (rnorm < 1e-6)
			goto GET_REACTION;

		if (rnorm != rnorm) {
			status = NO_INVERSE;
			goto EXIT;
		}

		status = solve_linear_system(K, perm, iperm, Kr, Lr, Ur,
					     residual, delta_disp);

		if (status != 0) {
			status = ELASTIC_SOLVER_FAILS;
			goto EXIT;
		}

		nb_vector_sum(2 * N_elems, displacement, delta_disp);

		assemble_global_damage(eval_dmg, glq, faces, D, rhs_damage);
		status = solve_linear_system(D, dmg_perm, dmg_iperm, Dr,
					     dmg_Lr, dmg_Ur, rhs_damage,
					     elem_damage);

		if (status != 0) {
			status = DAMAGE_SOLVER_FAILS;
			goto EXIT;
		}
		
		iter ++;
		if (iter >= max_iter) {
			status = MAX_ITER_REACHED;
			goto EXIT;
		}
		printf(" ----> DAMAGE ITER: %i (%e)\n", iter, rnorm);/* TEMP */
	}
GET_REACTION:
	nb_cvfa_assemble_global_stiffness(K, mesh, SMOOTH, intmsh, xc,
					  faces, material, analysis2D,
					  params2D, glq, eval_dmg);
	nb_sparse_multiply_vector(K, displacement, residual, 1);
	*monitor_reaction = nb_math_hypo(residual[id_elem_monitor * 2],
					 residual[id_elem_monitor*2+1]);
EXIT:
	if (SUCCESS_RESIDUAL_MIN != status)
		memcpy(displacement, aux_disp, N * sizeof(*aux_disp));

	printf(" >>>>> [%i] DAMAGE ITER: %i (%e) ... %e/%i\n",
	       status, iter, rnorm, bc_factor, max_iter);/* TEMPORAL */

	nb_bcond_finish(numeric_bcond);
	nb_sparse_destroy(Kr);
	nb_sparse_destroy(Dr);
	if (LU_ALLOCATION_FAILS != status) {
		nb_sparse_destroy(Lr);
		nb_sparse_destroy(Ur);
	}
	if (LU_ALLOCATION_FAILS_DMG != status) {
		nb_sparse_destroy(dmg_Lr);
		nb_sparse_destroy(dmg_Ur);
	}
	nb_soft_free_mem(bcond_size, numeric_bcond);
	return status;
}

static int solve_linear_system(const nb_sparse_t *A,
			       const uint32_t *perm,
			       const uint32_t *iperm,
			       nb_sparse_t *Ar,
			       nb_sparse_t *Lr,
			       nb_sparse_t *Ur,
			       double *b, double *x)
{
	uint32_t N = nb_sparse_get_size(A);
	nb_vector_permutation(N, b, perm, x);
	memcpy(b, x, N * sizeof(*b));

	nb_sparse_fill_permutation(A, Ar, perm, iperm);
	nb_sparse_decompose_LU(Ar, Lr, Ur, 2);
	nb_sparse_solve_LU(Lr, Ur, b, x);

	nb_vector_permutation(N, x, iperm, b);
	memcpy(x, b, N * sizeof(*x));
	return 0;
}

static void assemble_global_damage(const nb_cvfa_eval_damage_t * eval_dmg,
				   const nb_glquadrature_t *glq,
				   face_t **faces, nb_sparse_t *D, double *H)
{
	const eval_damage_data_t *dmg_data = eval_dmg->data;
	uint32_t N_elems = nb_mesh2D_get_N_elems(dmg_data->mesh);
	uint32_t N_faces = nb_mesh2D_get_N_edges(dmg_data->mesh);
	memset(H, 0, N_elems * sizeof(*H));
	nb_sparse_reset(D);
	
	if (MIDPOINT_VOL_INTEGRALS) {
		for (uint32_t i = 0; i < N_elems; i++)
			assemble_elem_damage(dmg_data, i, glq, D, H);
	}

	for (uint32_t i = 0; i < N_faces; i++)
		assemble_face_damage(dmg_data, faces[i], glq, D, H);
}

static void assemble_elem_damage(const eval_damage_data_t *dmg_data,
				 uint32_t elem_id,
				 const nb_glquadrature_t *glq,
				 nb_sparse_t *D, double *H)
{
	double area = nb_mesh2D_elem_get_area(dmg_data->mesh, elem_id);
	double energy = get_elem_energy(elem_id, dmg_data);
  	double h = nb_material_get_damage_length_scale(dmg_data->material);
	double G = nb_material_get_energy_release_rate(dmg_data->material);
	double val = 2 * energy * h/G;
	nb_sparse_add(D, elem_id, elem_id, area * (val + 1));
	H[elem_id] += area * val;
}

static void assemble_face_damage(const eval_damage_data_t *dmg_data,
				 face_t *face,
				 const nb_glquadrature_t *glq,
				 nb_sparse_t *D, double *H)
{
	if (nb_cvfa_face_is_internal(face, dmg_data->mesh)) {
		uint16_t N_sf = face->N_sf;
		for (uint16_t i = 0; i < N_sf; i++) {
			if (!MIDPOINT_VOL_INTEGRALS) {
				integrate_subvolume_damage(dmg_data, face,
							   i, glq, D, H);
			}
			integrate_subface_damage(dmg_data, face,
						 i, glq, D, H);
		}
	} else {
		integrate_exterior_subface_damage(dmg_data, face, glq, D, H);
	}
}

static void integrate_subvolume_damage(const eval_damage_data_t *dmg_data,
				       face_t *face, uint16_t subface_id,
				       const nb_glquadrature_t *glq,
				       nb_sparse_t *D, double *H)
{
	uint32_t id1 = face->elems[0];
	uint32_t id2 = face->elems[1];
	double c1[2], c2[2];
	c1[0] = dmg_data->xc[id1 * 2];
	c1[1] = dmg_data->xc[id1*2+1];
	c2[0] = dmg_data->xc[id2 * 2];
	c2[1] = dmg_data->xc[id2*2+1];

	const subface_t *subface = face->subfaces[subface_id];
	integrate_lateral_subvolume_damage(subface->x1, subface->x2, c1,
					   dmg_data, face, subface_id,
					   id1, glq, D, H);
	integrate_lateral_subvolume_damage(subface->x2, subface->x1, c2,
					   dmg_data, face, subface_id,
					   id2, glq, D, H);
}

static void integrate_lateral_subvolume_damage
				(const double st1[2],
				 const double st2[2],
				 const double st3[2],
				 const eval_damage_data_t *dmg_data,
				 const face_t *face, uint32_t subface_id,
				 uint32_t elem_id,
				 const nb_glquadrature_t *glq,
				 nb_sparse_t *D, double *H)
{
	double iJ[4];
	double detJ = nb_cvfa_subface_get_inverse_jacobian(st1, st2, st3, iJ);
	
	const subface_t *subface = face->subfaces[subface_id];
	for (uint16_t q1 = 0; q1 < glq->N; q1++) {
		double xaux[2];
		xaux[0] = glq->x[q1];
		for (uint16_t q2 = 0; q2 < glq->N; q2++) {
			xaux[1] = glq->x[q2];
			double sxi[2];
			sxi[0] = 0.5 * (1.0 + xaux[0]);
			sxi[1] = 0.25 * (1.0 - xaux[0]) * (1.0 + xaux[1]);
			double xq[2];
			nb_cvfa_get_interpolated_point(st1, st2, st3, sxi, xq);
			double wq = glq->w[q1] * glq->w[q2] *
				0.125 * (1.0 - xaux[0]) * detJ;
			if (nb_cvfa_subface_in_simplex(subface))
				integrate_svol_simplexwise_gp_dmg(dmg_data,
								  face,
								  subface_id,
								  elem_id,
								  wq, xq,
								  D, H);
			else
				integrate_svol_pairwise_gp_dmg(dmg_data,
							       face,
							       subface_id,
							       elem_id,
							       wq, xq, D, H);
		}
	}
}

static void integrate_svol_simplexwise_gp_dmg
				(const eval_damage_data_t *dmg_data,
				 const face_t *face, uint32_t subface_id,
				 uint32_t elem_id, double wq,
				 const double xq[2],
				 nb_sparse_t *D, double *H)
{
	const subface_t *subface = face->subfaces[subface_id];

	double t1[2], t2[2], t3[2];
	nb_cvfa_load_trg_points(dmg_data->intmsh,
				subface->trg_id, t1, t2, t3);

	double xi[2];
	nb_cvfa_get_normalized_point(dmg_data->smooth, t1, t2, t3, xq, xi);

	double energy = get_subface_energy(face, subface_id, xq, dmg_data);
  	double h = nb_material_get_damage_length_scale(dmg_data->material);
	double G = nb_material_get_energy_release_rate(dmg_data->material);
	double val = 2 * energy * h/G;
	for (uint8_t k = 0; k < 3; k++) {
		double Nk;
		if (k < 2)
			Nk = xi[k];
		else
			Nk = 1.0 - xi[0] - xi[1];

		uint32_t elem_k = nb_mesh2D_elem_get_adj(dmg_data->intmsh,
							 subface->trg_id, k);
		nb_sparse_add(D, elem_id, elem_k, wq * (val + 1) * Nk);
	}

	H[elem_id] += wq * val;
}

static void integrate_svol_pairwise_gp_dmg
				(const eval_damage_data_t *dmg_data,
				 const face_t *face, uint32_t subface_id,
				 uint32_t elem_id, double wq,
				 const double xq[2],
				 nb_sparse_t *D, double *H)
{
	uint32_t id1 = face->elems[0];
	uint32_t id2 = face->elems[1];
	const double *x1 = &(dmg_data->xc[id1 * 2]);
	const double *x2 = &(dmg_data->xc[id2 * 2]);
	double c[2];
	c[0] = x2[0] - x1[0];
	c[1] = x2[1] - x1[1];
	double cnorm2 = POW2(c[0]) + POW2(c[1]);
	c[0] /= cnorm2;
	c[1] /= cnorm2;

	double z = c[0] * (xq[0] - x1[0]) + c[1] * (xq[1] - x1[1]);

	double energy = get_subface_energy(face, subface_id, xq, dmg_data);
  	double h = nb_material_get_damage_length_scale(dmg_data->material);
	double G = nb_material_get_energy_release_rate(dmg_data->material);
	double val = 2 * energy * h/G;
	for (uint8_t k = 0; k < 2; k++) {
		double Nk;
		if (0 == k)
			Nk = 1.0 - z;
		else
			Nk = z;

		uint32_t elem_k = face->elems[k];
		nb_sparse_add(D, elem_id, elem_k, wq * (val + 1) * Nk);
	}

	H[elem_id] += wq * val;
}

static void integrate_subface_damage(const eval_damage_data_t *dmg_data,
				     face_t *face, uint16_t subface_id,
				     const nb_glquadrature_t *glq,
				     nb_sparse_t *D, double *H)
{
	subface_t *subface = face->subfaces[subface_id];
	for (uint8_t q = 0; q < glq->N; q++) {
		if (nb_cvfa_subface_in_simplex(subface))
			integrate_subface_simplexwise_gp_damage(dmg_data,
								face,
								subface_id,
								glq, q, D, H);
		else
			integrate_subface_pairwise_gp_damage(dmg_data,
							     face,
							     subface_id,
							     glq, q, D, H);
	}
}

static void integrate_subface_simplexwise_gp_damage
					(const eval_damage_data_t *dmg_data,
					 face_t *face, uint16_t subface_id,
					 const nb_glquadrature_t *glq,
					 uint8_t gp,
					 nb_sparse_t *D, double *H)
{
	subface_t *subface = face->subfaces[subface_id];
	double t1[2], t2[2], t3[2];
	nb_cvfa_load_trg_points(dmg_data->intmsh,
				subface->trg_id, t1, t2, t3);

	double xq[2];
	nb_cvfa_subface_get_xq(subface, glq, gp, xq);
	double xi[2];
	nb_cvfa_get_normalized_point(dmg_data->smooth, t1, t2, t3, xq, xi);

	double iJ[4];
	nb_cvfa_subface_get_inverse_jacobian(t1, t2, t3, iJ);

	double lf = nb_utils2D_get_dist(subface->x1, subface->x2);
	double wq = lf * glq->w[gp] * 0.5;

  	double h = nb_material_get_damage_length_scale(dmg_data->material);
	for (uint8_t k = 0; k < 3; k++) {
		double grad_xi[2];
		if (0 == k) {
			grad_xi[0] = 1;
			grad_xi[1] = 0;
		} else if (1 == k) {
			grad_xi[0] = 0;
			grad_xi[1] = 1;
		} else {
			grad_xi[0] = -1;
			grad_xi[1] = -1;
		}
		double grad[2];
		nb_cvfa_subface_get_grad(iJ, grad_xi, grad);

		double gradn = grad[0] * face->nf[0] + grad[1] * face->nf[1];

		uint32_t elem_k = nb_mesh2D_elem_get_adj(dmg_data->intmsh,
							 subface->trg_id, k);
		nb_sparse_add(D, face->elems[0], elem_k, -wq * POW2(h) * gradn);
		nb_sparse_add(D, face->elems[1], elem_k,  wq * POW2(h) * gradn);
	}
}

static void integrate_subface_pairwise_gp_damage
					(const eval_damage_data_t *dmg_data,
					 face_t *face, uint16_t subface_id,
					 const nb_glquadrature_t *glq,
					 uint8_t gp,
					 nb_sparse_t *D, double *H)
{
	uint32_t id1 = face->elems[0];
	uint32_t id2 = face->elems[1];
	const double *x1 = &(dmg_data->xc[id1 * 2]);
	const double *x2 = &(dmg_data->xc[id2 * 2]);
	double c[2];
	c[0] = x2[0] - x1[0];
	c[1] = x2[1] - x1[1];
	double cnorm2 = POW2(c[0]) + POW2(c[1]);
	c[0] /= cnorm2;
	c[1] /= cnorm2;

	const subface_t *subface = face->subfaces[subface_id];

	double xq[2];
	nb_cvfa_subface_get_xq(subface, glq, gp, xq);

	double lf = nb_utils2D_get_dist(subface->x1, subface->x2);
	double wq = lf * glq->w[gp] * 0.5;

  	double h = nb_material_get_damage_length_scale(dmg_data->material);
	for (uint8_t k = 0; k < 2; k++) {
		double dz;
		if (0 == k)
			dz = -1;
		else
			dz = 1;
		double grad[2];
		grad[0] = -dz * c[0] / cnorm2;
		grad[1] = -dz * c[1] / cnorm2;

		double gradn = grad[0] * face->nf[0] + grad[1] * face->nf[1];

		uint32_t elem_k = face->elems[k];
		nb_sparse_add(D, id1, elem_k, -wq * POW2(h) * gradn);
		nb_sparse_add(D, id2, elem_k,  wq * POW2(h) * gradn);
	}
}

static void integrate_exterior_subface_damage
					(const eval_damage_data_t *dmg_data,
					 face_t *face,
					 const nb_glquadrature_t *glq,
					 nb_sparse_t *D, double *H)
{
	/* TEMPORAL */
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
	case LU_ALLOCATION_FAILS:
		fprintf(stderr, "Allocation for LU decomposition fails");
		break;
	case LU_ALLOCATION_FAILS_DMG:
		fprintf(stderr, "Allocation for damage LU decomp. fails");
		break;
	default:
		fprintf(stderr, "Unkown error");		
	}
	fprintf(stderr, "\n");
}

static void get_reaction_log_name(const char *dir, char *reaction_log)
{
	sprintf(reaction_log, "%s/", dir);
	char *pch = strstr(reaction_log, "/");
	pch = pch + 1;
	sprintf(pch, "Reaction.log");
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

	fprintf(fp, "%i %e %e\n", iter, factor, reaction);
	fclose(fp);
EXIT:
	return;
}

void TEMPORAL_DRAW(const char *dir,
		   const nb_mesh2D_t *mesh, const double *elem_damage, 
		   double *face_damage, int i)
{
	uint32_t N_nodes = nb_mesh2D_get_N_nodes(mesh);
	double *damage = nb_allocate_mem(N_nodes * sizeof(*damage));

	nb_mesh2D_extrapolate_elems_to_nodes(mesh, 1, elem_damage, damage);
	char name[100];
	sprintf(name, "%s/CVFA_elem_dmg_%i.png", dir, i);
	nb_mesh2D_export_draw(mesh, name, 1000, 700, NB_NODE, NB_FIELD,
			      damage, true);
	nb_free_mem(damage);

	uint32_t N_faces = nb_mesh2D_get_N_edges(mesh);
	double max_dmg = 0;
	for (uint32_t i = 0; i < N_faces; i++) {
		if (max_dmg < face_damage[i])
			max_dmg = face_damage[i];
	}
	for (uint32_t i = 0; i < N_faces; i++)
		face_damage[i] /= max_dmg;
	nb_mesh2D_export_draw(mesh, "./CVFA_dmg.png", 1000, 800,
			      NB_FACE, NB_FIELD,
			      face_damage, true);
}

static void save_simulation(const char *dir,
			    const nb_mesh2D_t *mesh, const double *disp,
			    const double *elem_damage,
			    const double *face_damage, int iter, double time)
{
	char mesh_name[50];
	sprintf(mesh_name, "CVFA_DMG_mesh_%i.vtk", iter);
	char name[100];
	sprintf(name, "%s/%s", dir, mesh_name);
	
	nb_mesh2D_save_vtk(mesh, name);

	sprintf(name, "%s/CVFA_DMG_results_%i.log", dir, iter);
	save_data(name, mesh, disp, elem_damage, face_damage,
		  iter, time, mesh_name);

	TEMPORAL_DRAW(dir, mesh, elem_damage, (void*)face_damage, iter);	
}

static void save_data(const char *name,
		      const nb_mesh2D_t *mesh,
		      const double *disp,
		      const double *elem_damage,
		      const double *face_damage, int iter, double time,
		      const char *mesh_name)
{
	FILE *fp = fopen(name, "w");
	if (NULL == fp)
		goto EXIT;
	
	fprintf(fp, "# DAMAGE CVFA SIMULATION (cardoso.victore@gmail.com)\n");
	fprintf(fp, "# MESH_FILE\n %s\n", mesh_name);
	fprintf(fp, "# N_ITER\n %i\n", iter);
	fprintf(fp, "# TIME\n %e\n\n", time);

	uint32_t N_elems = nb_mesh2D_get_N_elems(mesh);
	fprintf(fp, "# N_ELEMS\n %i\n", N_elems);
	fprintf(fp, "# DISPLACEMENT\n");
	for (uint32_t i = 0; i < N_elems; i++)
		fprintf(fp, "%e %e\n", disp[i * 2], disp[i*2+1]);
	fprintf(fp, "# DAMAGE (ELEMENT INTEGRAL)\n");
	for (uint32_t i = 0; i < N_elems; i++)
		fprintf(fp, "%e\n", elem_damage[i]);

	uint32_t N_faces = nb_mesh2D_get_N_edges(mesh);
	fprintf(fp, "# N_FACES\n %i\n", N_faces);
	fprintf(fp, "# DAMAGE (FACE INTEGRAL)\n");
	for (uint32_t i = 0; i < N_faces; i++)
		fprintf(fp, "%e\n", face_damage[i]);

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

void nb_cvfa_compute_stress_from_damage_and_strain
					(const nb_mesh2D_t *mesh,
					 const nb_material_t *const material,
					 nb_analysis2D_t analysis2D,
					 const double* strain,
					 const double* damage,
					 double* stress /* Output */)
{
	uint32_t N_faces = nb_mesh2D_get_N_edges(mesh);
	for (uint32_t i = 0; i < N_faces; i++) {
		double D[4];
		nb_pde_get_constitutive_matrix(D, material, analysis2D);
		
		stress[i * 3] = (1 - damage[i]) *
			(strain[i * 3] * D[0] + strain[i*3+1] * D[1]);
		stress[i*3+1] = (1 - damage[i]) *
			(strain[i * 3] * D[1] + strain[i*3+1] * D[2]);
		stress[i*3+2] =  (1 - damage[i]) * strain[i*3+2] * D[3];
	}
}
