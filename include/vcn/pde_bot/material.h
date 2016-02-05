#ifndef __VCN_PDE_BOT_MATERIAL_H__
#define __VCN_PDE_BOT_MATERIAL_H__

typedef struct vcn_fem_material_s vcn_fem_material_t;

vcn_fem_material_t* vcn_fem_material_create();
void vcn_fem_material_destroy(vcn_fem_material_t* mat);

void vcn_fem_material_set_poisson_module(vcn_fem_material_t* mat, double pmodule);
void vcn_fem_material_set_elasticity_module(vcn_fem_material_t* mat, double emodule);
void vcn_fem_material_set_density(vcn_fem_material_t* mat, double density);
void vcn_fem_material_set_fracture_energy(vcn_fem_material_t* mat, double frac_energy);
void vcn_fem_material_set_traction_limit_stress(vcn_fem_material_t* mat, double max_stress);
void vcn_fem_material_set_compression_limit_stress(vcn_fem_material_t* mat, double max_stress);

double vcn_fem_material_get_poisson_module(const vcn_fem_material_t *const mat);
double vcn_fem_material_get_elasticity_module(const vcn_fem_material_t *const mat);
double vcn_fem_material_get_density(const vcn_fem_material_t *const mat);
double vcn_fem_material_get_fracture_energy(const vcn_fem_material_t *const mat);
double vcn_fem_material_get_traction_limit_stress(const vcn_fem_material_t *const mat);
double vcn_fem_material_get_compression_limit_stress(const vcn_fem_material_t *const mat);
void* vcn_fem_material_get_damage_function(const vcn_fem_material_t *const mat);

int vcn_fem_material_verify(vcn_fem_material_t* material);

#endif
