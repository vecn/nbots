#include <stdlib.h>
#include <stdbool.h>

#include "vcn/pde_bot/material.h"

struct vcn_fem_material_s {
	double poisson_module;
	double elasticity_module;
	double density;
	double fracture_energy;
	double traction_limit_stress;
	double compression_limit_stress;
	/* Damage function return the damage in
	 * the interval [0: No damage, 1: Full damage] */
	double (*damage)(const vcn_fem_material_t *const mat, 
			 double *strain,
			 double *r_damage_prev,
			 double *r_damage,
			 double characteristic_length_of_fractured_domain,
			 bool enable_plane_stress);
	double (*r0_damage)(const vcn_fem_material_t *const mat);
};

vcn_fem_material_t *vcn_fem_material_create(void)
{
	vcn_fem_material_t *mat = calloc(1, sizeof(*mat));
	return mat;
}

void vcn_fem_material_destroy(vcn_fem_material_t* mat){
	free(mat);
}

void vcn_fem_material_set_poisson_module(vcn_fem_material_t* mat, double pmodule){
	mat->poisson_module = pmodule;
}

void vcn_fem_material_set_elasticity_module(vcn_fem_material_t* mat, double emodule){
	mat->elasticity_module = emodule;
}

void vcn_fem_material_set_density(vcn_fem_material_t* mat, double density){
	mat->density = density;
}

void vcn_fem_material_set_fracture_energy(vcn_fem_material_t* mat, double frac_energy){
	mat->fracture_energy = frac_energy;
}

void vcn_fem_material_set_traction_limit_stress(vcn_fem_material_t* mat, double max_stress){
	mat->traction_limit_stress = max_stress;
}

void vcn_fem_material_set_compression_limit_stress(vcn_fem_material_t* mat, double max_stress){
	mat->compression_limit_stress = max_stress;
}

double vcn_fem_material_get_poisson_module(const vcn_fem_material_t *const mat){
	return mat->poisson_module;
}

double vcn_fem_material_get_elasticity_module(const vcn_fem_material_t *const mat){
	return mat->elasticity_module;
}

double vcn_fem_material_get_density(const vcn_fem_material_t *const mat){
	return mat->density;
}

double vcn_fem_material_get_fracture_energy(const vcn_fem_material_t *const mat){
	return mat->fracture_energy;
}

double vcn_fem_material_get_traction_limit_stress(const vcn_fem_material_t *const mat){
	return mat->traction_limit_stress;
}

double vcn_fem_material_get_compression_limit_stress(const vcn_fem_material_t *const mat){
	return mat->compression_limit_stress;
}

int vcn_fem_material_verify(vcn_fem_material_t* material){
	if(vcn_fem_material_get_poisson_module(material) < 0.0)
		return 1;
	if(vcn_fem_material_get_poisson_module(material) >= 0.5)
		return 2;
	return 0;
}
