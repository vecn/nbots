#ifndef __NB_PDE_BOT_CONTROL_VOLUME_ELASTICITY_2D_H__
#define __NB_PDE_BOT_CONTROL_VOLUME_ELASTICITY_2D_H__

#include <stdbool.h>
#include <stdint.h>

#include "nb/geometric_bot.h"
#include "nb/solver_bot.h"
#include "nb/pde_bot.h"

typedef struct {
	void *data;
	double (*get_damage)(const face_t *face, uint16_t subface_id,
			     uint8_t gp, const nb_glquadrature_t *glq,
			     const void *data);
} nb_cvfa_eval_damage_t;

void nb_cvfa_assemble_global_forces(double *F,
				    const nb_mesh2D_t *const mesh,
				    const nb_material_t *material,
				    bool enable_self_weight,
				    double gravity[2]);

void nb_cvfa_assemble_global_stiffness(nb_sparse_t *K,
				       const nb_mesh2D_t *const mesh,
				       int smooth,
				       const nb_mesh2D_t *intmsh,
				       const double *xc, face_t **faces,
				       const nb_material_t *material,
				       nb_analysis2D_t analysis2D,
				       nb_analysis2D_params *params2D,
				       const nb_glquadrature_t *glq,
				                     /* NULL for no damage */
				       const nb_cvfa_eval_damage_t* eval_dmg);

void nb_cvfa_get_normalized_point(int smooth, const double x1[2],
				  const double x2[2], const double x3[2],
				  const double xq[2], double xi[2]);

void nb_cvfa_get_interpolated_point(int smooth, const double x1[2],
				    const double x2[2], const double x3[2],
				    const double xi[2], double xq[2]);

void nb_cvfa_compute_strain(double *strain, char *boundary_mask,
			    face_t **faces, 
			    const nb_mesh2D_t *mesh, int smooth,
			    const nb_mesh2D_t *intmsh, const double *xc,
			    const nb_bcond_t *const bcond,
			    const double *disp,
			    const nb_glquadrature_t *glq);

void nb_cvfa_subface_get_strain(int smooth,
				const nb_mesh2D_t *intmsh,
				const face_t *face,
				const subface_t *subface,
				const double *xc,
				const double *disp,
				const nb_glquadrature_t *glq,
				uint8_t q, double *strain);

void nb_cvfa_subface_get_grad_strain(int smooth,
				     const nb_mesh2D_t *intmsh,
				     const face_t *face,
				     const subface_t *subface,
				     const double *xc,
				     const double *disp,
				     const nb_glquadrature_t *glq,
				     uint8_t q, double *grad_strain);

#endif
