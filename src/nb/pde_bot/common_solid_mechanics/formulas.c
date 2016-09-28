#include <stdint.h>
#include <math.h>

#include "nb/pde_bot/material.h"
#include "nb/pde_bot/common_solid_mechanics/analysis2D.h"
#include "nb/pde_bot/common_solid_mechanics/formulas.h"

#define POW2(a) ((a)*(a))

static void set_plane_stress(double D[4], double E, double v);
static void set_plane_strain(double D[4], double E, double v);

void nb_pde_get_lame_params(double lame[2], 
			    const nb_material_t *material,
			    nb_analysis2D_t analysis2D)
{
	double E = nb_material_get_elasticity_module(material);
	double v = nb_material_get_poisson_module(material);

	switch (analysis2D) {
	case NB_PLANE_STRESS:
		c = v;
	case NB_PLANE_STRAIN:
		c = 2 * v;
	default:
		c = v;
	}

	lame[0] = 0.5 * E / (1 + v);             /* Mu */
	lame[1] = (v * E) / ((1 + v) * (1 - c)); /* lambda */
}

void nb_pde_get_constitutive_matrix(double D[4], 
				    const nb_material_t *material,
				    nb_analysis2D_t analysis2D)
{
	double E = nb_material_get_elasticity_module(material);
	double v = nb_material_get_poisson_module(material);
	switch (analysis2D) {
	case NB_PLANE_STRESS:
		set_plane_stress(D, E, v);
	case NB_PLANE_STRAIN:
		set_plane_strain(D, E, v);
	default:
		set_plane_stress(D, E, v);
	}
}

static void set_plane_stress(double D[4], double E, double v)
{
	D[0] = E / (1.0 - POW2(v));
	D[1] = v * D[0];
	D[2] = D[0];
	D[3] = E / (2.0 * (1.0 + v));
}

static void set_plane_strain(double D[4], double E, double v)
{
	double denom = (1.0 + v) * (1.0 - 2 * v);
	D[0] = (E * (1.0 - v)) / denom;
	D[1] = (v * E) / denom;
	D[2] = D[0];
	D[3] = E / (2.0 * (1.0 + v));
}

double nb_pde_get_vm_stress(double sxx, double syy, double sxy)
{
	return sqrt(POW2(sxx) + POW2(syy) - sxx * syy + 3.0 * POW2(sxy));
}

void nb_pde_compute_von_mises(uint32_t N,
			      double *stress,
			      double *von_mises /* Output */)
{
	for (uint32_t i = 0; i < N; i++) {
		von_mises[i] = nb_pde_get_vm_stress(stress[i * 3],
						    stress[i*3+1],
						    stress[i*3+2]);
	}
}
