#include <stdlib.h>
#include <stdbool.h>

#include "nb/memory_bot.h"
#include "nb/pde_bot/material.h"

struct nb_material_s {
	double elasticity_module;/* TEMPORAL: Modulus instead of module */
	double poisson_module;
	double density;
	double fracture_energy;
	double traction_limit_stress;
	double compression_limit_stress;
	double energy_release_rate;
	double damage_length_scale;
	/* Damage function return the damage in
	 * the interval [0: No damage, 1: Full damage] */
	double (*damage)(const nb_material_t *const mat, 
			 double *strain,
			 double *r_damage_prev,
			 double *r_damage,
			 double characteristic_length_of_fractured_domain,
			 bool enable_plane_stress);
	double (*r0_damage)(const nb_material_t *const mat);
};

nb_material_t *nb_material_create(void)
{
	nb_material_t *mat = nb_allocate_zero_mem(sizeof(*mat));
	return mat;
}

nb_material_t *nb_material_clone(nb_material_t* mat)
{
	nb_material_t *clone = nb_material_create();
	clone->elasticity_module = mat->elasticity_module;
	clone->poisson_module = mat->poisson_module;
	clone->density = mat->density;
	clone->fracture_energy = mat->fracture_energy;
	clone->traction_limit_stress = mat->traction_limit_stress;
	clone->compression_limit_stress = mat->compression_limit_stress;
	clone->damage = mat->damage;
	clone->r0_damage = mat->r0_damage;
	return clone;
}

void nb_material_destroy(nb_material_t* mat)
{
	nb_free_mem(mat);
}

void nb_material_set_poisson_module(nb_material_t* mat, double pmodule)
{
	mat->poisson_module = pmodule;
}

void nb_material_set_elasticity_module(nb_material_t* mat, double emodule)
{
	mat->elasticity_module = emodule;
}

void nb_material_set_density(nb_material_t* mat, double density)
{
	mat->density = density;
}

void nb_material_set_fracture_energy(nb_material_t* mat, double frac_energy)
{
	mat->fracture_energy = frac_energy;
}

void nb_material_set_traction_limit_stress(nb_material_t* mat,
					   double max_stress)
{
	mat->traction_limit_stress = max_stress;
}

void nb_material_set_compression_limit_stress(nb_material_t* mat,
					      double max_stress)
{
	mat->compression_limit_stress = max_stress;
}

void nb_material_set_energy_release_rate(nb_material_t *mat, double energy_rr)
{
	mat->energy_release_rate = energy_rr;
}

void nb_material_set_damage_length_scale(nb_material_t *mat,
					 double length_scale)
{
	mat->damage_length_scale = length_scale;
}

double nb_material_get_poisson_module(const nb_material_t *const mat)
{
	return mat->poisson_module;
}

double nb_material_get_elasticity_module(const nb_material_t *const mat)
{
	return mat->elasticity_module;
}

double nb_material_get_density(const nb_material_t *const mat)
{
	return mat->density;
}

double nb_material_get_fracture_energy(const nb_material_t *const mat)
{
	return mat->fracture_energy;
}

double nb_material_get_traction_limit_stress(const nb_material_t *const mat)
{
	return mat->traction_limit_stress;
}

double nb_material_get_compression_limit_stress(const nb_material_t *const mat)
{
	return mat->compression_limit_stress;
}

double nb_material_get_energy_release_rate(const nb_material_t *const mat)
{
	return mat->energy_release_rate;
}

double nb_material_get_damage_length_scale(const nb_material_t *const mat)
{
	return mat->damage_length_scale;
}

int nb_material_verify(nb_material_t* material)
{
	if(nb_material_get_poisson_module(material) < 0.0)
		return 1;
	if(nb_material_get_poisson_module(material) >= 0.5)
		return 2;
	return 0;
}
