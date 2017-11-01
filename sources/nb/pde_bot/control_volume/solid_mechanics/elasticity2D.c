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
static void get_jacobian(const double t1[2],
			 const double t2[2],
			 const double t3[2],
			 double J[4]);
static void subface_get_normalized_grad(int smooth, uint8_t i,
					const double xi[2],
					double grad_xi[2]);
static double get_spline(int smooth, double x);
static double get_deriv_spline(int smooth, double x);
static double get_second_deriv_spline(int smooth, double x);
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
static void subface_get_strain_simplexwise(int smooth,
					   const nb_mesh2D_t *intmsh,
					   const subface_t *subface,
					   const double *disp,
					   const double xq[2], double *strain);
static void subface_get_strain_pairwise
			(int smooth,
			 const face_t *face,
			 const subface_t *subface,
			 const double *xc,
			 const double *disp,
			 const double xq[2],
			 double *strain);
static void get_boundary_face_strain(face_t **faces, uint32_t face_id,
				     const nb_bcond_t *bcond,
				     const double *disp, double *strain);

static void subface_get_grad_strain_simplexwise
			(int smooth,
			 const nb_mesh2D_t *intmsh,
			 const subface_t *subface,
			 const double *disp,
			 const double xq[2],
			 double *grad_strain);
static void subface_get_grad_strain_pairwise
			(int smooth,
			 const face_t *face,
			 const subface_t *subface,
			 const double *xc,
			 const double *disp,
			 const double xq[2],
			 double *grad_strain);

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
	if (nb_cvfa_subface_in_simplex(subface)) {
		for (uint8_t q = 0; q < glq->N; q++)
			integrate_subface_simplexwise(K, mesh, smooth, intmsh,
						      face, subface_id, D,
						      params2D, glq, q,
						      eval_dmg);
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
			nb_vector_scale(12, Kf, POW2(1.0 - damage));

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
	nb_cvfa_subface_get_xq(subface, glq, q, xq);

	double xi[2];
	nb_cvfa_get_normalized_point(t1, t2, t3, xq, xi);
	
	double wq = lf * glq->w[q] * 0.5;

	double iJ[4];
	nb_cvfa_subface_get_inverse_jacobian(t1, t2, t3, iJ);

	double factor = wq * params2D->thickness; /* TEMP Check PLANE_STRAIN */
	for (uint8_t i = 0; i < 3; i++) {
		double grad_xi[2];
		subface_get_normalized_grad(smooth, i, xi, grad_xi);
		double grad[2];
		nb_cvfa_subface_get_grad(iJ, grad_xi, grad);
		double Kfi[4];
		subface_get_nodal_contribution(D, face->nf, grad, Kfi);
		Kf[i * 2] = factor * Kfi[0];
		Kf[i*2+1] = factor * Kfi[1];
		Kf[6 + i * 2] = factor * Kfi[2];
		Kf[6 + i*2+1] = factor * Kfi[3];
	}
}

double nb_cvfa_subface_get_inverse_jacobian(const double t1[2],
					    const double t2[2],
					    const double t3[2],
					    double iJ[4])
{
	get_jacobian(t1, t2, t3, iJ);

	double aux = iJ[1]; /* Transpose Jacobian */
	iJ[1] = iJ[2];
	iJ[2] = aux;

	double det = nb_matrix_2X2_inverse_destructive(iJ);

	return det;
}

static void get_jacobian(const double t1[2],
			 const double t2[2],
			 const double t3[2],
			 double J[4])
{
	/* Jacobian = D_{psi} x*/
	J[0] = t1[0] - t3[0];
	J[1] = t2[0] - t3[0];
	J[2] = t1[1] - t3[1];
	J[3] = t2[1] - t3[1];
}

static void subface_get_normalized_grad(int smooth, uint8_t i,
					const double xi[2],
					double grad_xi[2])
{
	if (0 == i) {
		double dPx = get_deriv_spline(smooth, xi[0]);
		grad_xi[0] = dPx;
		grad_xi[1] = 0;
	} else if (1 == i) {
		double dPy = get_deriv_spline(smooth, xi[1]);
		grad_xi[0] = 0;
		grad_xi[1] = dPy;
	} else {
		double dPx = get_deriv_spline(smooth, xi[0]);
		double dPy = get_deriv_spline(smooth, xi[1]);
		grad_xi[0] = -dPx;
		grad_xi[1] = -dPy;
	}
}

void nb_cvfa_get_normalized_point(const double x1[2],
				  const double x2[2], const double x3[2],
				  const double xq[2], double xi[2])
{
	double Jd[4];
	get_jacobian(x1, x2, x3, Jd);
	
	double b[2];
	b[0] = xq[0] - x3[0];
	b[1] = xq[1] - x3[1];

	nb_matrix_2X2_inverse_destructive(Jd);

	xi[0] = Jd[0] * b[0] + Jd[1] * b[1];
	xi[1] = Jd[2] * b[0] + Jd[3] * b[1];
}

void nb_cvfa_get_interpolated_point(const double x1[2],
				    const double x2[2], const double x3[2],
				    const double xi[2], double xq[2])
{
	double P1 = xi[0];
	double P2 = xi[1];
	double P3 = 1.0 - P1 - P2;
	xq[0] = x1[0] * P1 + x2[0] * P2 + x3[0] * P3;
	xq[1] = x1[1] * P1 + x2[1] * P2 + x3[1] * P3;
}

void nb_cvfa_get_interpolated_disp(int smooth, const double u1[2],
				   const double u2[2], const double u3[2],
				   const double xi[2], double uq[2])
{
	double P1 = get_spline(smooth, xi[0]);
	double P2 = get_spline(smooth, xi[1]);
	double P3 = 1.0 - P1 - P2;
	uq[0] = u1[0] * P1 + u2[0] * P2 + u3[0] * P3;
	uq[1] = u1[1] * P1 + u2[1] * P2 + u3[1] * P3;
}

static double get_spline(int smooth, double x)
{
	double spline;
	void* aux = &spline;
	switch (smooth) {
	case 0:
		spline = x;
		break;
	case 1:
		spline = 3 * POW2(x) - 2 * POW3(x);
		break;
	case 2:
		spline = 10 * POW3(x) - 15 * pow(x, 4) + 6 * pow(x, 5);
		break;
	case 3:
		spline = 35 * pow(x, 4) - 84 * pow(x, 5) +
			70 * pow(x, 6) - 20 * pow(x, 7);
		break;
	case 4:
		spline = 126 * pow(x, 5) - 420 * pow(x, 6) +
			540 * pow(x, 7) - 315 * pow(x, 8) + 70 * pow(x, 9);
		break;
	case 5:
		spline = 462 * pow(x, 6) - 1980 * pow(x, 7) + 3465 * pow(x, 8) -
			3080 * pow(x, 9) + 1386 * pow(x, 10) - 252 * pow(x, 11);
		break;
	default:
		spline = 1716 * pow(x, 7) - 9009 * pow(x, 8) +
			20020 * pow(x, 9) - 24024 * pow(x, 10) +
			16380 * pow(x, 11) - 6006 * pow(x, 12) +
			924 * pow(x, 13);
	}
	memcpy(aux, &x, 8);/* Plan HTUMZ */
	return spline;
}

static double get_deriv_spline(int smooth, double x)
{
	double deriv;
	void* aux = &deriv;
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
	*((double*)aux) = 1;/* Plan HTUMZ */
	return deriv;
}

static double get_second_deriv_spline(int smooth, double x)
{
	double deriv;
	void* aux = &deriv;
	switch (smooth) {
	case 0:
		deriv = 0;
		break;
	case 1:
		deriv = 6 - 12 * x;
		break;
	case 2:
		deriv = 60 * x - 180 * POW2(x) + 120 * POW3(x);
		break;
	case 3:
		deriv = 420 * POW2(x) - 1680 * pow(x, 3) +
			2100 * pow(x, 4) - 840 * pow(x, 5);
		break;
	case 4:
		deriv = 2520 * pow(x, 3) - 12600 * pow(x, 4) +
			22680 * pow(x, 5) - 17640 * pow(x, 6) +
			5040 * pow(x, 7);
		break;
	case 5:
		deriv = 13860 * pow(x, 4) - 83160 * pow(x, 5) +
			194040 * pow(x, 6) - 221760 * pow(x, 7) +
			124740 * pow(x, 8) - 27720 * pow(x, 9);
		break;
	default:
		/* Smooth = 6 */
		deriv = 72072 * pow(x, 5) - 504504 * pow(x, 6) +
			1441440 * pow(x, 7) - 2162160 * pow (x, 8) +
			1801800 * pow(x, 9) - 792792 * pow(x, 10) +
			144144  * pow(x, 11);
	}
	memset(aux, 0, 8);/* Plan HTUMZ */
	return deriv;

}

void nb_cvfa_subface_get_grad(const double iJ[4], const double grad_xi[2],
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
			nb_vector_scale(8, Kf, POW2(1.0 - damage));

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
	nb_cvfa_subface_get_xq(subface, glq, q, xq);
	
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
			double xq[2];
			nb_cvfa_subface_get_xq(subface, glq, q, xq);
			double gp_strain[3];
			nb_cvfa_subface_get_strain(smooth, intmsh, face,
						   subface, xc, disp, xq,
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

void nb_cvfa_subface_get_xq(const subface_t *subface,
			    const nb_glquadrature_t *glq,
			    uint8_t q, double xq[2])
{
	double xstep = (glq->x[q] + 1) / 2.0;
	xq[0] = subface->x1[0] + xstep * (subface->x2[0] - subface->x1[0]);
	xq[1] = subface->x1[1] + xstep * (subface->x2[1] - subface->x1[1]);
}

void nb_cvfa_subface_get_strain(int smooth,
				const nb_mesh2D_t *intmsh,
				const face_t *face,
				const subface_t *subface,
				const double *xc,
				const double *disp,
				const double xq[2],
				double *strain)
{
	if (nb_cvfa_subface_in_simplex(subface))
		subface_get_strain_simplexwise(smooth, intmsh, subface,
					       disp, xq, strain);
	else
		subface_get_strain_pairwise(smooth, face, subface, xc,
					    disp, xq, strain);
}

static void subface_get_strain_simplexwise(int smooth,
					   const nb_mesh2D_t *intmsh,
					   const subface_t *subface,
					   const double *disp,
					   const double xq[2], double *strain)
{
	double t1[2], t2[2], t3[2];
	nb_cvfa_load_trg_points(intmsh, subface->trg_id, t1, t2, t3);

	double xi[2];
	nb_cvfa_get_normalized_point(t1, t2, t3, xq, xi);

	double iJ[4];
	nb_cvfa_subface_get_inverse_jacobian(t1, t2, t3, iJ);

	memset(strain, 0, 3 * sizeof(*strain));
	for (uint8_t i = 0; i < 3; i++) {
		double grad_xi[2];
		subface_get_normalized_grad(smooth, i, xi, grad_xi);
		double grad[2];
		nb_cvfa_subface_get_grad(iJ, grad_xi, grad);

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
			 const double xq[2],
			 double *strain)
{
	uint32_t id1 = face->elems[0];
	uint32_t id2 = face->elems[1];
	double c1[2], c2[2];
	c1[0] = xc[id1 * 2];
	c1[1] = xc[id1*2+1];
	c2[0] = xc[id2 * 2];
	c2[1] = xc[id2*2+1];
	
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

void nb_cvfa_subface_get_grad_strain(int smooth,
				     const nb_mesh2D_t *intmsh,
				     const face_t *face,
				     const subface_t *subface,
				     const double *xc,
				     const double *disp,
				     const double xq[2],
				     double *grad_strain)
{
	if (nb_cvfa_subface_in_simplex(subface))
		subface_get_grad_strain_simplexwise(smooth, intmsh, subface,
						    disp, xq, grad_strain);
	else
		subface_get_grad_strain_pairwise(smooth, face, subface, xc,
						 disp, xq, grad_strain);
}

static void subface_get_grad_strain_simplexwise
			(int smooth,
			 const nb_mesh2D_t *intmsh,
			 const subface_t *subface,
			 const double *disp,
			 const double xq[2],
			 double *grad_strain)
{
	double t1[2], t2[2], t3[2];
	nb_cvfa_load_trg_points(intmsh, subface->trg_id, t1, t2, t3);

	double xi[2];
	nb_cvfa_get_normalized_point(t1, t2, t3, xq, xi);

	uint32_t e1 = nb_mesh2D_elem_get_adj(intmsh, subface->trg_id, 0);
	double u1 = disp[e1 * 2];
	double v1 = disp[e1 * 2 + 1];
	uint32_t e2 = nb_mesh2D_elem_get_adj(intmsh, subface->trg_id, 1);
	double u2 = disp[e2 * 2];
	double v2 = disp[e2 * 2 + 1];
	uint32_t e3 = nb_mesh2D_elem_get_adj(intmsh, subface->trg_id, 2);
	double u3 = disp[e3 * 2];
	double v3 = disp[e3 * 2 + 1];

	double d2hx = get_second_deriv_spline(smooth, xi[0]);
	double d2hy = get_second_deriv_spline(smooth, xi[1]);

	grad_strain[0] = d2hx * (u1 - u3);
	grad_strain[1] = 0;
	grad_strain[2] = 0;
	grad_strain[3] = d2hy * (v2 - v3);
	grad_strain[4] = d2hx * (v1 - v3);
	grad_strain[5] = d2hy * (u2 - u3);
}

static void subface_get_grad_strain_pairwise
			(int smooth,
			 const face_t *face,
			 const subface_t *subface,
			 const double *xc,
			 const double *disp,
			 const double xq[2],
			 double *grad_strain)
{
	uint32_t id1 = face->elems[0];
	uint32_t id2 = face->elems[1];
	double c[2];
	c[0] = xc[id2 * 2] - xc[id1 * 2];
	c[1] = xc[id2*2+1] - xc[id1*2+1];
	double cnorm2 = POW2(c[0]) + POW2(c[1]);
	c[0] /= cnorm2;
	c[1] /= cnorm2;

	double z = c[0] * (xq[0] - xc[id1 * 2]) + c[1] * (xq[1] - xc[id1*2+1]);
	double d2hz = get_second_deriv_spline(smooth, z);

	double u1 = disp[id1 * 2];
	double u2 = disp[id1*2+1];
	double v1 = disp[id2 * 2];
	double v2 = disp[id2*2+1];

	grad_strain[0] = d2hz * POW2(c[0]) * (u2 - u1);
	grad_strain[1] = d2hz * c[0] *c[1] * (u2 - u1);
	grad_strain[2] = d2hz * c[0] *c[1] * (v2 - v1);
	grad_strain[3] = d2hz * POW2(c[1]) * (v2 - v1);
	grad_strain[4] = d2hz * c[0] *c[1] * (u2 - u1) +
		         d2hz * POW2(c[0]) * (v2 - v1);
	grad_strain[5] = d2hz * POW2(c[1]) * (u2 - u1) +
		         d2hz * c[0] *c[1] * (v2 - v1);
}
