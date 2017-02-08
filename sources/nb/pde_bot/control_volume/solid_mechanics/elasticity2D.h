#ifndef __NB_PDE_BOT_CONTROL_VOLUME_ELASTICITY_2D_H__
#define __NB_PDE_BOT_CONTROL_VOLUME_ELASTICITY_2D_H__

#include <stdbool.h>
#include <stdint.h>

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

double nb_cvfa_subface_get_inverse_jacobian(const double t1[2],
					    const double t2[2],
					    const double t3[2],
					    double iJ[4]);

void nb_cvfa_subface_get_grad(const double iJ[4], const double grad_xi[2],
			      double grad[2]);
void nb_cvfa_get_normalized_point(const double x1[2],
				  const double x2[2], const double x3[2],
				  const double xq[2], double xi[2]);

void nb_cvfa_get_interpolated_point(const double x1[2],
				    const double x2[2], const double x3[2],
				    const double xi[2], double xq[2]);

void nb_cvfa_get_interpolated_disp(int smooth, const double u1[2],
				   const double u2[2], const double u3[2],
				   const double xi[2], double uq[2]);

void nb_cvfa_compute_strain(double *strain, char *boundary_mask,
			    face_t **faces, 
			    const nb_mesh2D_t *mesh, int smooth,
			    const nb_mesh2D_t *intmsh, const double *xc,
			    const nb_bcond_t *const bcond,
			    const double *disp,
			    const nb_glquadrature_t *glq);

void nb_cvfa_subface_get_xq(const subface_t *subface,
			    const nb_glquadrature_t *glq,
			    uint8_t q, double xq[2]);

void nb_cvfa_subface_get_strain(int smooth,
				const nb_mesh2D_t *intmsh,
				const face_t *face,
				const subface_t *subface,
				const double *xc,
				const double *disp,
				const double xq[2],
				double *strain);

void nb_cvfa_subface_get_grad_strain(int smooth,
				     const nb_mesh2D_t *intmsh,
				     const face_t *face,
				     const subface_t *subface,
				     const double *xc,
				     const double *disp,
				     const double xq[2],
				     double *grad_strain);

#endif
