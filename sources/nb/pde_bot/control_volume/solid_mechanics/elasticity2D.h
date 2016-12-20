#ifndef __NB_PDE_BOT_CONTROL_VOLUME_ELASTICITY2D_H__
#define __NB_PDE_BOT_CONTROL_VOLUME_ELASTICITY2D_H__

#include <stdbool.h>

#include "nb/geometric_bot.h"
#include "nb/pde_bot.h"

typedef struct subface_s subface_t;

typedef struct {
	double nf[2];
	uint32_t elems[2];
	double x1[2], x2[2];
	uint16_t N_sf;
	subface_t **subfaces;
} face_t;

struct subface_s {
	uint8_t N_int;/*     Zero: Pairwise      */
                      /* Not zero: Simplex-wise  */
	double x1[2], x2[2];
	uint32_t trg_id;
};

int nb_cvfa_solve_elasticity_equation(const nb_mesh2D_t *mesh,
				      const nb_material_t *material,
				      const nb_bcond_t *bcond,
				      int smooth,
				      bool enable_self_weight,
				      double gravity[2],
				      nb_analysis2D_t analysis2D,
				      nb_analysis2D_params *params2D,
				      double *displacement, /* Output */
				      const nb_mesh2D_t *intmsh,
				      const double *xc,
				      face_t **faces,
				      double *F, nb_sparse_t *K,
				      nb_glquadrature_t *glq);

void nb_cvfa_compute_strain(double *strain, char *boundary_mask,
			    face_t **faces,
			    const nb_mesh2D_t *const mesh, int smooth,
			    const nb_mesh2D_t *intmsh, const double *xc,
			    const nb_bcond_t *const bcond,
			    const double *disp,
			    const nb_glquadrature_t *glq);

#endif
