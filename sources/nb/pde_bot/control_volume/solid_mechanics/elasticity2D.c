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

#include "../integration_mesh.h"

#include "elasticity2D.h"
#include "set_bconditions.h"

#define INV3 0.33333333333333333333334
#define POW2(a) ((a)*(a))
#define POW3(a) ((a)*(a)*(a))


static void integrate_elem_force(const nb_mesh2D_t *mesh,
				 const nb_material_t *material,
				 bool enable_self_weight,
				 double gravity[2],
				 uint32_t elem_id,
				 double *F);
static void assemble_face(nb_sparse_t *K,
			  const nb_mesh2D_t *const mesh,
			  int smooth,
			  const nb_mesh2D_t *intmsh,
			  const double *xc, face_t *face,
			  const nb_material_t *material,
			  nb_analysis2D_t analysis2D,
			  nb_analysis2D_params *params2D,
			  const nb_glquadrature_t *glq,
			                /* NULL for no damage */
			  const nb_cvfa_eval_damage_t* eval_dmg);
static void integrate_subface(nb_sparse_t *K,
			      const nb_mesh2D_t *const mesh,
			      int smooth,
			      const nb_mesh2D_t *intmsh,
			      const double *xc, face_t *face,
			      const double D[4],
			      nb_analysis2D_params *params2D,
			      uint16_t subface_id,
			      const nb_glquadrature_t *glq,
			                    /* NULL for no damage */
			      const nb_cvfa_eval_damage_t* eval_dmg);
static void integrate_subface_simplexwise(nb_sparse_t *K,
					  const nb_mesh2D_t *const mesh,
					  int smooth,
					  const nb_mesh2D_t *intmsh,
					  face_t *face, uint16_t subface_id,
					  const double D[4],
					  nb_analysis2D_params *params2D,
					  const nb_glquadrature_t *glq,
					  uint8_t q,
					                /* NULL for no damage */
					  const nb_cvfa_eval_damage_t* eval_dmg);
static void integrate_Kf(const nb_mesh2D_t *const mesh, int smooth,
			 const nb_mesh2D_t *intmsh, face_t *face,
			 uint16_t subface_id, const double D[4],
			 nb_analysis2D_params *params2D, double Kf[12],
			 const nb_glquadrature_t *glq, uint8_t q);
static double subface_get_inverse_jacobian(int smooth,
					   const double t1[2],
					   const double t2[2],
					   const double t3[2],
					   double iJ[4],
					   const double xi[2]);
static void get_jacobian(int somooth,
			 const double t1[2],
			 const double t2[2],
			 const double t3[2],
			 double J[4],
			 const double xi[2]);
static void subface_get_normalized_grad(int smooth, uint8_t i,
					const double xi[2],
					double grad_xi[2]);
static double get_deriv_spline(int smooth, double x);
static double get_spline_inv(int smooth, double x);
static void subface_get_grad(const double iJ[4], const double grad_xi[2],
			     double grad[2]);
static void subface_get_nodal_contribution(const double D[4],
					   const double nf[2],
					   const double grad[2],
					   double Kfi[4]);
static void add_Kf_to_K(face_t *face, const nb_mesh2D_t *intmsh,
			uint16_t subface_id, const double Kf[12],
			nb_sparse_t *K);
static void integrate_subface_pairwise(nb_sparse_t *K,
				       const nb_mesh2D_t *const mesh,
				       int smooth,
				       const double *xc, face_t *face,
				       uint16_t subface_id,
				       const double D[4],
				       nb_analysis2D_params *params2D,
				       const nb_glquadrature_t *glq,
				       uint8_t q,
				                     /* NULL for no damage */
				       const nb_cvfa_eval_damage_t* eval_dmg);
static void integrate_Kf_pairwise(const nb_mesh2D_t *const mesh, int smooth,
				  const double *xc, face_t *face,
				  uint16_t subface_id, const double D[4],
				  nb_analysis2D_params *params2D, double Kf[8],
				  const nb_glquadrature_t *glq, uint8_t q);
static void face_get_grad_pairwise(int smooth,
				   const double c1[2], const double c2[2],
				   double grad[2], const double x[2]);
static void add_Kf_to_K_pairwise(face_t *face, const double Kf[8],
				 nb_sparse_t *K);
static void get_face_strain(face_t **faces, uint32_t face_id,
			    const nb_mesh2D_t *mesh, int smooth,
			    const nb_mesh2D_t *intmsh,
			    const double *xc,
			    const nb_bcond_t *const bcond,
			    const double *disp,
			    double *strain,
			    char *boundary_mask,
			    const nb_glquadrature_t *glq);
static void get_internal_face_strain(face_t **faces, uint32_t face_id,
				     int smooth,
				     const nb_mesh2D_t *intmsh,
				     const double *xc,
				     const double *disp, double *strain,
				     const nb_glquadrature_t *glq);
static void subface_get_strain_simplexwise
			(int smooth,
			 const nb_mesh2D_t *intmsh,
			 const subface_t *subface,
			 const double *disp,
			 const nb_glquadrature_t *glq,
			 uint8_t q, double *strain);
static void subface_get_strain_pairwise
			(int smooth,
			 const face_t *face,
			 const subface_t *subface,
			 const double *xc,
			 const double *disp,
			 const nb_glquadrature_t *glq,
			 uint8_t q, double *strain);
static void get_boundary_face_strain(face_t **faces, uint32_t face_id,
				     const nb_bcond_t *bcond,
				     const double *disp, double *strain);

void nb_cvfa_assemble_global_forces(double *F,
				    const nb_mesh2D_t *const mesh,
				    const nb_material_t *material,
				    bool enable_self_weight,
				    double gravity[2])
{
	uint32_t N_elems = nb_mesh2D_get_N_elems(mesh);
	memset(F, 0, N_elems * 2 * sizeof(*F));
	for (uint32_t i = 0; i < N_elems; i++) {
		integrate_elem_force(mesh, material, enable_self_weight,
				     gravity, i, F);
	}
}

static void integrate_elem_force(const nb_mesh2D_t *mesh,
				 const nb_material_t *material,
				 bool enable_self_weight,
				 double gravity[2],
				 uint32_t elem_id,
				 double *F)
{
	if (enable_self_weight) {
		double area = nb_mesh2D_elem_get_area(mesh, elem_id);
		double mass = area * nb_material_get_density(material);
		F[elem_id * 2] += mass * gravity[0];
		F[elem_id*2+1] += mass * gravity[1];
	}
}

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
				       const nb_cvfa_eval_damage_t* eval_dmg)
{
	nb_sparse_reset(K);
	uint32_t N_faces = nb_mesh2D_get_N_edges(mesh);
	for (uint32_t i = 0; i < N_faces; i++) {
		assemble_face(K, mesh, smooth, intmsh, xc, faces[i], material,
			      analysis2D, params2D, glq, eval_dmg);
	}
}

static void assemble_face(nb_sparse_t *K,
			  const nb_mesh2D_t *const mesh,
			  int smooth,
			  const nb_mesh2D_t *intmsh,
			  const double *xc, face_t *face,
			  const nb_material_t *material,
			  nb_analysis2D_t analysis2D,
			  nb_analysis2D_params *params2D,
			  const nb_glquadrature_t *glq,
			               /* NULL for no damage */
			  const nb_cvfa_eval_damage_t* eval_dmg)
{	
	double D[4];
	nb_pde_get_constitutive_matrix(D, material, analysis2D);

	if (nb_cvfa_face_is_internal(face, mesh)) {
		uint16_t N_sf = face->N_sf;
		for (uint16_t i = 0; i < N_sf; i++) {
			integrate_subface(K, mesh, smooth, intmsh, xc, face,
					  D, params2D, i, glq, eval_dmg);
		}
	}
}

static void integrate_subface(nb_sparse_t *K,
			      const nb_mesh2D_t *const mesh,
			      int smooth,
			      const nb_mesh2D_t *intmsh,
			      const double *xc, face_t *face,
			      const double D[4],
			      nb_analysis2D_params *params2D,
			      uint16_t subface_id,
			      const nb_glquadrature_t *glq,
			                    /* NULL for no damage */
			      const nb_cvfa_eval_damage_t* eval_dmg)
{
	subface_t *subface = face->subfaces[subface_id];
	if (subface->N_int > 0) {
		for (uint8_t q = 0; q < glq->N; q++)
			integrate_subface_simplexwise(K, mesh, smooth, intmsh,
						      face, subface_id, D,
						      params2D, glq, q, eval_dmg);
	} else {
		for (uint8_t q = 0; q < glq->N; q++)
			integrate_subface_pairwise(K, mesh, smooth, xc, face,
						   subface_id, D, params2D,
						   glq, q, eval_dmg);
	}
}

static void integrate_subface_simplexwise(nb_sparse_t *K,
					  const nb_mesh2D_t *const mesh,
					  int smooth,
					  const nb_mesh2D_t *intmsh,
					  face_t *face, uint16_t subface_id,
					  const double D[4],
					  nb_analysis2D_params *params2D,
					  const nb_glquadrature_t *glq,
					  uint8_t q,
					                /* NULL for no damage */
					  const nb_cvfa_eval_damage_t* eval_dmg)
{
	double damage = 0.0;
	if (NULL != eval_dmg)
		damage = eval_dmg->get_damage(face, subface_id,
					      q, glq, eval_dmg->data);
	if (damage < 1.0) {
		double Kf[12];
		integrate_Kf(mesh, smooth, intmsh, face, subface_id, D,
			     params2D, Kf, glq, q);

		if (damage > 0.0)
			nb_vector_scale(12, Kf, 1.0 - damage);

		add_Kf_to_K(face, intmsh, subface_id, Kf, K);
	}
}

static void integrate_Kf(const nb_mesh2D_t *const mesh, int smooth,
			 const nb_mesh2D_t *intmsh, face_t *face,
			 uint16_t subface_id, const double D[4],
			 nb_analysis2D_params *params2D, double Kf[12],
			 const nb_glquadrature_t *glq,
			 uint8_t q)
{
	subface_t *subface = face->subfaces[subface_id];

	double t1[2], t2[2], t3[2];
	nb_cvfa_load_trg_points(intmsh, subface->trg_id, t1, t2, t3);

	double lf = nb_utils2D_get_dist(subface->x1, subface->x2);

	double xq[2];
	double xstep = (glq->x[q] + 1) / 2.0;
	xq[0] = subface->x1[0] + xstep * (subface->x2[0] - subface->x1[0]);
	xq[1] = subface->x1[1] + xstep * (subface->x2[1] - subface->x1[1]);
	double xi[2];
	nb_cvfa_get_normalized_point(smooth, t1, t2, t3, xq, xi);
	
	double wq = lf * glq->w[q] * 0.5;

	double iJ[4];
	subface_get_inverse_jacobian(smooth, t1, t2, t3, iJ, xi);

	double factor = wq * params2D->thickness;
	for (uint8_t i = 0; i < 3; i++) {
		double grad_xi[2];
		subface_get_normalized_grad(smooth, i, xi, grad_xi);
		double grad[2];
		subface_get_grad(iJ, grad_xi, grad);
		double Kfi[4];
		subface_get_nodal_contribution(D, face->nf, grad, Kfi);
		Kf[i * 2] = factor * Kfi[0];
		Kf[i*2+1] = factor * Kfi[1];
		Kf[6 + i * 2] = factor * Kfi[2];
		Kf[6 + i*2+1] = factor * Kfi[3];
	}
}

static double subface_get_inverse_jacobian(int smooth,
					   const double t1[2],
					   const double t2[2],
					   const double t3[2],
					   double iJ[4],
					   const double xi[2])
{
	get_jacobian(smooth, t1, t2, t3, iJ, xi);

	double aux = iJ[1];
	iJ[1] = iJ[2];
	iJ[2] = aux;

	double det = nb_matrix_2X2_inverse_destructive(iJ);

	return det;
}

static void get_jacobian(int smooth,
			 const double t1[2],
			 const double t2[2],
			 const double t3[2],
			 double J[4],
			 const double xi[2])
{
	/* Jacobian = D_{psi} x*/
	if (0 == smooth) {
		J[0] = t2[0] - t1[0];
		J[1] = t3[0] - t1[0];
		J[2] = t2[1] - t1[1];
		J[3] = t3[1] - t1[1];
	} else {
		memset(J, 0, 4 * sizeof(*J));
		double grad_xi[2];
		subface_get_normalized_grad(smooth, 0, xi, grad_xi);
		J[0] += grad_xi[0] * t1[0];
		J[1] += grad_xi[1] * t1[0];
		J[2] += grad_xi[0] * t1[1];
		J[3] += grad_xi[1] * t1[1];

		subface_get_normalized_grad(smooth, 1, xi, grad_xi);
		J[0] += grad_xi[0] * t2[0];
		J[1] += grad_xi[1] * t2[0];
		J[2] += grad_xi[0] * t2[1];
		J[3] += grad_xi[1] * t2[1];

		subface_get_normalized_grad(smooth, 2, xi, grad_xi);
		J[0] += grad_xi[0] * t3[0];
		J[1] += grad_xi[1] * t3[0];
		J[2] += grad_xi[0] * t3[1];
		J[3] += grad_xi[1] * t3[1];
	}
}

static void subface_get_normalized_grad(int smooth, uint8_t i,
					const double xi[2],
					double grad_xi[2])
{
	if (0 == i) {
		double dPx = get_deriv_spline(smooth, xi[0]);
		double dPy = get_deriv_spline(smooth, xi[1]);
		grad_xi[0] = -dPx;
		grad_xi[1] = -dPy;
	} else if (1 == i) {
		double dPx = get_deriv_spline(smooth, xi[0]);
		grad_xi[0] = dPx;
		grad_xi[1] = 0;
	} else {
		double dPy = get_deriv_spline(smooth, xi[1]);
		grad_xi[0] = 0;
		grad_xi[1] = dPy;
	}
}

void nb_cvfa_get_normalized_point(int smooth, const double x1[2],
				  const double x2[2], const double x3[2],
				  const double xq[2], double xi[2])
{
	double Jd[4];
	Jd[0] = x2[0] - x1[0];
	Jd[1] = x3[0] - x1[0];
	Jd[2] = x2[1] - x1[1];
	Jd[3] = x3[1] - x1[1];
	
	double b[2];
	b[0] = xq[0] - x1[0];
	b[1] = xq[1] - x1[1];

	nb_matrix_2X2_inverse_destructive(Jd);

	xi[0] = get_spline_inv(smooth, Jd[0] * b[0] + Jd[1] * b[1]);
	xi[1] = get_spline_inv(smooth, Jd[2] * b[0] + Jd[3] * b[1]);
}

static double get_deriv_spline(int smooth, double x)
{
	double deriv;
	switch (smooth) {
	case 0:
		deriv = 1;
		break;
	case 1:
		deriv = 6 * x - 6 * POW2(x);
		break;
	case 2:
		deriv = 30 * POW2(x) - 60 * POW3(x) + 30 * pow(x, 4);
		break;
	case 3:
		deriv = 140 * pow(x, 3) - 420 * pow(x, 4) +
			420 * pow(x, 5) - 140 * pow(x, 6);
		break;
	case 4:
		deriv = 630 * pow(x, 4) - 2520 * pow(x, 5) + 3780 * pow(x, 6) -
			2520 * pow(x, 7) + 630 * pow(x, 8);
		break;
	case 5:
		deriv = 2772 * pow(x, 5) - 13860 * pow(x, 6) +
			27720 * pow(x, 7) - 27720 * pow(x, 8) +
			13860 * pow(x, 9) - 2772 * pow(x, 10);
		break;
	default:
		/* Smooth = 6 */
		deriv = 12012 * pow(x, 6) - 72072 * pow(x, 7) +
			180180 * pow(x, 8) - 240240 * pow (x, 9) +
			180180 * pow(x, 10) - 72072 * pow(x, 11) +
			12012  * pow(x, 12);
	}
	return deriv;
}

static double get_spline_inv(int smooth, double x)
{
	double inv;
	switch (smooth) {
	case 0:
		inv = x;
		break;
	case 1:
		inv = 1.66667 * x - 2 * POW2(x) + 1.33334 * POW3(x);
		break;
	case 2:
		inv = 1.9333 * x - 2.8 * POW2(x) + (28.0/15.0) * POW3(x);
		break;
	case 3:
		inv = 2.0867 * x - 3.2571 * POW2(x) + 2.1714 * POW3(x);
		break;
	case 4:
		inv = 2.1873 * x - 3.5619 * POW2(x) + 2.3746 * POW3(x);
		break;
	case 5:
		inv = 2.26 * x - 3.7825 * POW2(x) + 2.5223 * POW3(x);
		break;
	default:
		inv = 2.3180 * x - 7.908 * POW2(x) + 7.908 * POW3(x);
	}
	return inv;
}


static void subface_get_grad(const double iJ[4], const double grad_xi[2],
			     double grad[2])
{
	grad[0] = iJ[0] * grad_xi[0] + iJ[1] * grad_xi[1];
	grad[1] = iJ[2] * grad_xi[0] + iJ[3] * grad_xi[1];
}

static void subface_get_nodal_contribution(const double D[4],
					   const double nf[2],
					   const double grad[2],
					   double Kfi[4])
{
	Kfi[0] = grad[0] * nf[0] * D[0] + grad[1] * nf[1] * D[3];
	Kfi[1] = grad[1] * nf[0] * D[1] + grad[0] * nf[1] * D[3];
	Kfi[2] = grad[0] * nf[1] * D[1] + grad[1] * nf[0] * D[3];
	Kfi[3] = grad[1] * nf[1] * D[2] + grad[0] * nf[0] * D[3];
}

static void add_Kf_to_K(face_t *face, const nb_mesh2D_t *intmsh,
			uint16_t subface_id, const double Kf[12],
			nb_sparse_t *K)
{
	uint32_t i = face->elems[0];
	uint32_t j = face->elems[1];
	uint32_t trg_id = face->subfaces[subface_id]->trg_id;
	for (uint8_t m = 0; m < 3; m++) {
		uint32_t k = nb_mesh2D_elem_get_adj(intmsh, trg_id, m);
		nb_sparse_add(K, i * 2, k * 2, -Kf[m * 2]);
		nb_sparse_add(K, i * 2, k*2+1, -Kf[m*2+1]);
		nb_sparse_add(K, i*2+1, k * 2, -Kf[6 + m * 2]);
		nb_sparse_add(K, i*2+1, k*2+1, -Kf[6 + m*2+1]);

		nb_sparse_add(K, j * 2, k * 2, Kf[m * 2]);
		nb_sparse_add(K, j * 2, k*2+1, Kf[m*2+1]);
		nb_sparse_add(K, j*2+1, k * 2, Kf[6 + m * 2]);
		nb_sparse_add(K, j*2+1, k*2+1, Kf[6 + m*2+1]);
	}
}

static void integrate_subface_pairwise(nb_sparse_t *K,
				       const nb_mesh2D_t *const mesh,
				       int smooth,
				       const double *xc, face_t *face,
				       uint16_t subface_id,
				       const double D[4],
				       nb_analysis2D_params *params2D,
				       const nb_glquadrature_t *glq,
				       uint8_t q,
				                     /* NULL for no damage */
				       const nb_cvfa_eval_damage_t* eval_dmg)
{
	double damage = 0.0;
	if (NULL != eval_dmg)
		damage = eval_dmg->get_damage(face, subface_id,
					      q, glq, eval_dmg->data);
	if (damage < 1.0) {
		double Kf[8];
		integrate_Kf_pairwise(mesh, smooth, xc, face, subface_id, D,
				      params2D, Kf, glq, q);
		if (damage > 0.0)
			nb_vector_scale(8, Kf, 1.0 - damage);

		add_Kf_to_K_pairwise(face, Kf, K);
	}
}

static void integrate_Kf_pairwise(const nb_mesh2D_t *const mesh, int smooth,
				  const double *xc, face_t *face,
				  uint16_t subface_id, const double D[4],
				  nb_analysis2D_params *params2D, double Kf[8],
				  const nb_glquadrature_t *glq, uint8_t q)
{
	uint32_t id1 = face->elems[0];
	uint32_t id2 = face->elems[1];
	double c1[2], c2[2];
	c1[0] = xc[id1 * 2];
	c1[1] = xc[id1*2+1];
	c2[0] = xc[id2 * 2];
	c2[1] = xc[id2*2+1];

	subface_t *subface = face->subfaces[subface_id];
	double lf = nb_utils2D_get_dist(subface->x1, subface->x2);

	double xq[2];
	double xstep = (glq->x[q] + 1) / 2.0;
	xq[0] = subface->x1[0] + xstep * (subface->x2[0] - subface->x1[0]);
	xq[1] = subface->x1[1] + xstep * (subface->x2[1] - subface->x1[1]);
	
	double wq = lf * glq->w[q] * 0.5;

	double factor = wq * params2D->thickness;
	for (uint8_t i = 0; i < 2; i++) {
		double grad[2];
		if (0 == i)
			face_get_grad_pairwise(smooth, c1, c2, grad, xq);
		else
			face_get_grad_pairwise(smooth, c2, c1, grad, xq);
		double Kfi[4];
		subface_get_nodal_contribution(D, face->nf, grad, Kfi);
		Kf[i * 2] = factor * Kfi[0];
		Kf[i*2+1] = factor * Kfi[1];
		Kf[4 + i * 2] = factor * Kfi[2];
		Kf[4 + i*2+1] = factor * Kfi[3];
	}
}

static void face_get_grad_pairwise(int smooth,
				   const double c1[2], const double c2[2],
				   double grad[2], const double x[2])
{
	double xdiff = c2[0] - c1[0];
	double ydiff = c2[1] - c1[1];
	double d2 = nb_utils2D_get_dist2(c1, c2);

	double dot = (x[0] - c1[0]) * xdiff + (x[1] - c1[1]) * ydiff;
	double z = dot / d2;

	double dPz = get_deriv_spline(smooth, z);

	grad[0] = -dPz * xdiff / d2;
	grad[1] = -dPz * ydiff / d2;
}

static void add_Kf_to_K_pairwise(face_t *face, const double Kf[8],
				 nb_sparse_t *K)
{
	uint32_t i = face->elems[0];
	uint32_t j = face->elems[1];
	for (uint8_t m = 0; m < 2; m++) {
		uint32_t k = face->elems[m];
		nb_sparse_add(K, i * 2, k * 2, -Kf[m * 2]);
		nb_sparse_add(K, i * 2, k*2+1, -Kf[m*2+1]);
		nb_sparse_add(K, i*2+1, k * 2, -Kf[4 + m * 2]);
		nb_sparse_add(K, i*2+1, k*2+1, -Kf[4 + m*2+1]);

		nb_sparse_add(K, j * 2, k * 2, Kf[m * 2]);
		nb_sparse_add(K, j * 2, k*2+1, Kf[m*2+1]);
		nb_sparse_add(K, j*2+1, k * 2, Kf[4 + m * 2]);
		nb_sparse_add(K, j*2+1, k*2+1, Kf[4 + m*2+1]);
	}
}

void nb_cvfa_compute_strain(double *strain, char *boundary_mask,
			    face_t **faces,
			    const nb_mesh2D_t *mesh, int smooth,
			    const nb_mesh2D_t *intmsh, const double *xc,
			    const nb_bcond_t *const bcond,
			    const double *disp,
			    const nb_glquadrature_t *glq)
{
	uint32_t N_faces = nb_mesh2D_get_N_edges(mesh);

 	for (uint32_t i = 0; i < N_faces; i++) {
		get_face_strain(faces, i, mesh, smooth, intmsh, xc, bcond,
				disp, strain, boundary_mask, glq);
	}
}

static void get_face_strain(face_t **faces, uint32_t face_id,
			    const nb_mesh2D_t *mesh, int smooth,
			    const nb_mesh2D_t *intmsh,
			    const double *xc,
			    const nb_bcond_t *const bcond,
			    const double *disp,
			    double *strain,
			    char *boundary_mask,
			    const nb_glquadrature_t *glq)
{
	if (nb_cvfa_face_is_internal(faces[face_id], mesh)) {
		boundary_mask[face_id] = 0;
		get_internal_face_strain(faces, face_id, smooth,
					 intmsh, xc, disp, strain, glq);
	} else {
		boundary_mask[face_id] = 1;
		get_boundary_face_strain(faces, face_id, bcond, disp, strain);
	}
}

static void get_internal_face_strain(face_t **faces, uint32_t face_id,
				     int smooth,
				     const nb_mesh2D_t *intmsh,
				     const double *xc,
				     const double *disp, double *strain,
				     const nb_glquadrature_t *glq)
{
	
	memset(&(strain[face_id*3]), 0, 3 * sizeof(*strain));
	face_t *face = faces[face_id];
	for (uint16_t i = 0; i < face->N_sf; i++) {
		subface_t *subface = face->subfaces[i];
		double lf = nb_utils2D_get_dist(subface->x1, subface->x2);
		for (uint8_t q = 0; q < glq->N; q++) {
			double gp_strain[3];
			nb_cvfa_subface_get_strain(smooth, intmsh, face,
						   subface, xc, disp, glq, q,
						   gp_strain);
			double wq = lf * glq->w[q] * 0.5;
			double factor = wq /* * params2D->thickness*/;/* TEMP */
			strain[face_id * 3] += factor * gp_strain[0];
			strain[face_id*3+1] += factor * gp_strain[1];
			strain[face_id*3+2] += factor * gp_strain[2];
		}
	}
	double length = nb_utils2D_get_dist(face->x1, face->x2);
	strain[face_id * 3] /= length;
	strain[face_id*3+1] /= length;
	strain[face_id*3+2] /= length;
}

void nb_cvfa_subface_get_strain(int smooth,
				const nb_mesh2D_t *intmsh,
				const face_t *face,
				const subface_t *subface,
				const double *xc,
				const double *disp,
				const nb_glquadrature_t *glq,
				uint8_t q, double *strain)
{
	if (subface->N_int > 0)
		subface_get_strain_simplexwise(smooth, intmsh, subface,
					       disp, glq, q, strain);
	else
		subface_get_strain_pairwise(smooth, face, subface, xc,
					    disp, glq, q, strain);
}

static void subface_get_strain_simplexwise
			(int smooth,
			 const nb_mesh2D_t *intmsh,
			 const subface_t *subface,
			 const double *disp,
			 const nb_glquadrature_t *glq,
			 uint8_t q, double *strain)
{
	double t1[2], t2[2], t3[2];
	nb_cvfa_load_trg_points(intmsh, subface->trg_id, t1, t2, t3);

	double xq[2];
	double xstep = (glq->x[q] + 1) / 2.0;
	xq[0] = subface->x1[0] + xstep * (subface->x2[0] - subface->x1[0]);
	xq[1] = subface->x1[1] + xstep * (subface->x2[1] - subface->x1[1]);
	double xi[2];
	nb_cvfa_get_normalized_point(smooth, t1, t2, t3, xq, xi);

	double iJ[4];
	subface_get_inverse_jacobian(smooth, t1, t2, t3, iJ, xi);

	memset(strain, 0, 3 * sizeof(*strain));
	for (uint8_t i = 0; i < 3; i++) {
		double grad_xi[2];
		subface_get_normalized_grad(smooth, i, xi, grad_xi);
		double grad[2];
		subface_get_grad(iJ, grad_xi, grad);

		uint32_t elem_id =
			nb_mesh2D_elem_get_adj(intmsh, subface->trg_id, i);

		double u = disp[elem_id * 2];
		double v = disp[elem_id * 2 + 1];

		strain[0] += (grad[0] * u);
		strain[1] += (grad[1] * v);
		strain[2] += (grad[1] * u + grad[0] * v);
	}
}

static void subface_get_strain_pairwise
			(int smooth,
			 const face_t *face,
			 const subface_t *subface,
			 const double *xc,
			 const double *disp,
			 const nb_glquadrature_t *glq,
			 uint8_t q, double *strain)
{
	uint32_t id1 = face->elems[0];
	uint32_t id2 = face->elems[1];
	double c1[2], c2[2];
	c1[0] = xc[id1 * 2];
	c1[1] = xc[id1*2+1];
	c2[0] = xc[id2 * 2];
	c2[1] = xc[id2*2+1];

	double xq[2];
	double xstep = (glq->x[q] + 1) / 2.0;
	xq[0] = subface->x1[0] + xstep * (subface->x2[0] - subface->x1[0]);
	xq[1] = subface->x1[1] + xstep * (subface->x2[1] - subface->x1[1]);
	
	memset(strain, 0, 3 * sizeof(*strain));
	for (uint8_t i = 0; i < 2; i++) {
		double grad[2];
		if (0 == i)
			face_get_grad_pairwise(smooth, c1, c2, grad, xq);
		else
			face_get_grad_pairwise(smooth, c2, c1, grad, xq);

		uint32_t elem_id = face->elems[i];
		double u = disp[elem_id * 2];
		double v = disp[elem_id * 2 + 1];

		strain[0] += (grad[0] * u);
		strain[1] += (grad[1] * v);
		strain[2] += (grad[1] * u + grad[0] * v);
	}
}

static void get_boundary_face_strain(face_t **faces, uint32_t face_id,
				     const nb_bcond_t *bcond,
				     const double *disp, double *strain)
{
	memset(&(strain[face_id * 3]), 0, 3 * sizeof(*strain));
}
