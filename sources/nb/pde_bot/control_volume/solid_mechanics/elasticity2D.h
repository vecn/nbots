#ifndef __NB_PDE_BOT_CONTROL_VOLUME_ELASTICITY_2D_H__
#define __NB_PDE_BOT_CONTROL_VOLUME_ELASTICITY_2D_H__

#include <stdbool.h>

#include "nb/geometric_bot.h"
#include "nb/solver_bot.h"
#include "nb/pde_bot.h"

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
				       const nb_glquadrature_t *glq);

void nb_cvfa_compute_strain(double *strain, char *boundary_mask,
			    face_t **faces,
			    const nb_mesh2D_t *const mesh, int smooth,
			    const nb_mesh2D_t *intmsh, const double *xc,
			    const nb_bcond_t *const bcond,
			    const double *disp,
			    const nb_glquadrature_t *glq);

#endif
