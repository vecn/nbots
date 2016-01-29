/******************************************************************************
 *   FEM Bot: Finit element Method Bot                                        *
 *   2011-2015 Victor Eduardo Cardoso Nungaray                                *
 *   Twitter: @victore_cardoso                                                *
 *   email: victorc@cimat.mx                                                  *
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "vcn/math_bot.h"
#include "vcn/cfreader_cat.h"
#include "vcn/eigen_bot.h"
#include "vcn/container_bot.h"
#include "vcn/graph_bot.h"
#include "vcn/fem_bot.h"

/************************ Structure definition ****************************/
struct vcn_fem_elem_s{
	uint32_t N_nodes;
	uint32_t N_Gauss_points;
	double *psi, *eta;          /* Normalized space coordinates of Gauss points */
	double * gp_weight;         /* Integration weights of the Gauss points */
	double (**Ni)(double psi, double eta);       /* Shape functions */
	double (**dNi_dpsi)(double psi, double eta); /* Shape functions derivatives */
	double (**dNi_deta)(double psi, double eta); /* Shape functions derivatives */
};

struct vcn_fem_material_s{
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

struct vcn_fem_implicit_s{
	uint32_t N_steps;
	uint32_t N_max_iter;
	uint32_t N_max_iter_without_enhance;
	double residual_tolerance;
};

/**************** Shape functions and derivatives *************************/
/*          Type: Triangular element */
/* Interpolation: Linear             */
double N1(double psi, double eta){return 1.0 - psi - eta;}
double N2(double psi, double eta){return psi;}
double N3(double psi, double eta){return eta;}
double dN1_dpsi(double psi, double eta){return -1.0;}
double dN1_deta(double psi, double eta){return -1.0;}
double dN2_dpsi(double psi, double eta){return 1.0;}
double dN2_deta(double psi, double eta){return 0.0;}
double dN3_dpsi(double psi, double eta){return 0.0;}
double dN3_deta(double psi, double eta){return 1.0;}

/************************ Public functions definition **********************/
vcn_fem_implicit_t* vcn_fem_implicit_create(){
	return (vcn_fem_implicit_t*)
		calloc(1, sizeof(vcn_fem_implicit_t));
}

void vcn_fem_implicit_destroy(vcn_fem_implicit_t* isparams){
	free(isparams);
}

void vcn_fem_implicit_set_N_steps(vcn_fem_implicit_t* isparams,
				  uint32_t N_steps){
	isparams->N_steps = N_steps;
}

void vcn_fem_implicit_set_N_max_iter(vcn_fem_implicit_t* isparams,
				     uint32_t N_max_iter){
	isparams->N_max_iter = N_max_iter;
}

void vcn_fem_implicit_set_N_max_iter_without_enhance
(vcn_fem_implicit_t* isparams, uint32_t N_max_iter){
	isparams->N_max_iter_without_enhance = N_max_iter;
}

void vcn_fem_implicit_set_residual_tolerance
(vcn_fem_implicit_t* isparams, double penergy_tol)
{
	isparams->residual_tolerance = penergy_tol;
}

uint32_t vcn_fem_implicit_get_N_steps(vcn_fem_implicit_t* isparams){
	return isparams->N_steps;
}

uint32_t vcn_fem_implicit_get_N_max_iter(vcn_fem_implicit_t* isparams){
	return isparams->N_max_iter;
}

uint32_t vcn_fem_implicit_get_N_max_iter_without_enhance
(vcn_fem_implicit_t* isparams){
	return isparams->N_max_iter_without_enhance;
}

double vcn_fem_implicit_get_residual_tolerance
(vcn_fem_implicit_t* isparams)
{
	return isparams->residual_tolerance;
}

double tension_damage_r0(const vcn_fem_material_t *const mat){
	return vcn_fem_material_get_traction_limit_stress(mat) /
		sqrt(vcn_fem_material_get_elasticity_module(mat));
}

double tension_damage
(const vcn_fem_material_t *const mat, 
 double *strain,
 double *r_damage_prev,
 double *r_damage,
 double characteristic_length_of_fractured_domain,
 bool enable_plane_stress)
{
	/* Material properties */
	double E = vcn_fem_material_get_elasticity_module(mat);
	double v = vcn_fem_material_get_poisson_module(mat);
	double Gf = vcn_fem_material_get_fracture_energy(mat);
	double ft = vcn_fem_material_get_traction_limit_stress(mat);

	/* Get constitutive matrix */
	double d11 = E/(1.0 - vcn_math_pow2(v));
	double d12 = v*d11;
    
	if(!enable_plane_stress){
		d11 = (E*(1.0-v))/((1.0 + v)*(1.0-2*v));
		d12 = (v*d11)/(1.0-v);
	}    
	double d22 = d11;
	double d33 = E/(2.0*(1.0+v));

	/* Compute effective stress */
	double effective_stress[3];
	effective_stress[0] = d11 * strain[0] + d12 * strain[1];
	effective_stress[1] = d12 * strain[0] + d22 * strain[1];
	effective_stress[2] = d33 * strain[2];

	/* Compute principal stress using Mohr's circle */
	double sigma_avg = 0.5 * (effective_stress[0] + effective_stress[1]);
	double R = sqrt(0.25 * 
			vcn_math_pow2(effective_stress[0] - effective_stress[1]) +
			vcn_math_pow2(effective_stress[2]));
	double positive_stress1 = vcn_math_maxd(0, sigma_avg + R);
	double positive_stress2 = vcn_math_maxd(0, sigma_avg - R);

	/* Compute inverse of D */
	double detD2x2 = d11*d22 - d12*d12;
	double id11 =  d22/detD2x2;
	double id12 = -d12/detD2x2;
	double id22 =  d11/detD2x2;
	/* Compute stress ^T strain */
	double sTs =
		id11 * vcn_math_pow2(positive_stress1) + 
		2 * id12 * positive_stress1 * positive_stress2 +
		id22 * vcn_math_pow2(positive_stress2);

	/* Compute tau */
	double tau = sqrt(sTs);

	/* Compute and return damage */
	double r0 = mat->r0_damage(mat);
	r_damage[0] = vcn_math_maxd(r_damage_prev[0], tau);
	double div =  
		(Gf/characteristic_length_of_fractured_domain)*(E/vcn_math_pow2(ft));
	double A = 1.0 / (div - 0.5);
	double G = 1.0 - (r0/r_damage[0])*exp(A*(1.0-(r_damage[0]/r0)));
	return G;
}

double tension_truncated_damage_r0(const vcn_fem_material_t *const mat){
	return vcn_fem_material_get_traction_limit_stress(mat);
}

double tension_truncated_damage
(const vcn_fem_material_t *const mat, 
 double *strain,
 double *r_damage_prev,
 double* r_damage,
 double characteristic_length_of_fractured_domain,
 bool enable_plane_stress)
{
	/* Material properties */
	double E = vcn_fem_material_get_elasticity_module(mat);
	double v = vcn_fem_material_get_poisson_module(mat);
	double Gf = vcn_fem_material_get_fracture_energy(mat);
	double ft = vcn_fem_material_get_traction_limit_stress(mat);

	/* Get constitutive matrix */
	double d11 = E/(1.0 - vcn_math_pow2(v));
	double d12 = v*d11;
    
	if(!enable_plane_stress){
		d11 = (E*(1.0-v))/((1.0 + v)*(1.0-2*v));
		d12 = (v*d11)/(1.0-v);
	}    
	double d22 = d11;
	double d33 = E/(2.0*(1.0+v));

	/* Compute effective stress */
	double effective_stress[3];
	effective_stress[0] = d11 * strain[0] + d12 * strain[1];
	effective_stress[1] = d12 * strain[0] + d22 * strain[1];
	effective_stress[2] = d33 * strain[2] * 0.5;

	/* Compute principal stress using Mohr's circle */
	double sigma_avg = 0.5 * (effective_stress[0] + effective_stress[1]);
	double R = sqrt(0.25 * 
			vcn_math_pow2(effective_stress[0] - effective_stress[1]) +
			vcn_math_pow2(effective_stress[2]));
	double main_stress1 = sigma_avg + R;
	double main_stress2 = sigma_avg - R;

	/* Compute tau */
	double tau = vcn_math_maxd(0, main_stress1) + 
		vcn_math_maxd(0, main_stress2);

	/* Compute and return damage */
	double r0 = mat->r0_damage(mat);
	double beta = 2.0;
	r_damage[0] = vcn_math_maxd(r_damage_prev[0], tau);
	double pi = (beta * characteristic_length_of_fractured_domain *
		     vcn_math_pow2(ft))
		/(2*E*Gf);
	if(pi > 1.0) pi = 1.0;
	double Hs = pi/(1.0 - pi);
	double G = 1.0 - (r0/r_damage[0])*exp(2.0*Hs*(1.0-r_damage[0]/r0));
	return G;
}

vcn_fem_material_t* vcn_fem_vcn_fem_material_create(){
	vcn_fem_material_t* mat = (vcn_fem_material_t*)calloc(1, sizeof(vcn_fem_material_t));
	mat->damage = tension_truncated_damage; /* OPPORTUNITY: Dynamic assignment */
	mat->r0_damage = tension_truncated_damage_r0;
	return mat;
}

void vcn_fem_vcn_fem_material_destroy(vcn_fem_material_t* mat){
	free(mat);
}

void vcn_fem_vcn_fem_material_set_poisson_module(vcn_fem_material_t* mat, double pmodule){
	mat->poisson_module = pmodule;
}

void vcn_fem_vcn_fem_material_set_elasticity_module(vcn_fem_material_t* mat, double emodule){
	mat->elasticity_module = emodule;
}

void vcn_fem_vcn_fem_material_set_density(vcn_fem_material_t* mat, double density){
	mat->density = density;
}

void vcn_fem_vcn_fem_material_set_fracture_energy(vcn_fem_material_t* mat, double frac_energy){
	mat->fracture_energy = frac_energy;
}

void vcn_fem_vcn_fem_material_set_traction_limit_stress(vcn_fem_material_t* mat, double max_stress){
	mat->traction_limit_stress = max_stress;
}

void vcn_fem_vcn_fem_material_set_compression_limit_stress(vcn_fem_material_t* mat, double max_stress){
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

void* vcn_fem_material_get_damage_function(const vcn_fem_material_t *const mat){
	return mat->damage;
}

int vcn_fem_material_verify(vcn_fem_material_t* material){
	if(vcn_fem_material_get_poisson_module(material) < 0.0)
		return 1;
	if(vcn_fem_material_get_poisson_module(material) >= 0.5)
		return 2;
	return 0;
}

vcn_bcond_t* vcn_fem_bcond_create(){
	vcn_bcond_t* bconditions =
		(vcn_bcond_t*)calloc(1, sizeof(vcn_bcond_t));
	return bconditions;
}

vcn_bcond_t* vcn_fem_bcond_clone(const vcn_bcond_t*const bconditions){
	vcn_bcond_t* clone =
		(vcn_bcond_t*)calloc(1, sizeof(vcn_bcond_t));
	vcn_fem_bcond_copy(clone, bconditions);
	return clone;
}

vcn_bcond_t* vcn_fem_bcond_read(const char* filename){
	vcn_bcond_t* bconditions =
		(vcn_bcond_t*)calloc(1, sizeof(vcn_bcond_t));
  
	/* Initialize custom format to read file */
	vcn_cfreader_t* cfreader = vcn_cfreader_create(filename, "#");
	if(cfreader == NULL){
		vcn_fem_bcond_destroy(bconditions);
		return NULL;
	}
	/* Read number of Dirichlet conditions upon vertices */
	if(vcn_cfreader_read_uint32_t(cfreader, &(bconditions->N_Dirichlet_on_vtx)) != 0){
		vcn_cfreader_destroy(cfreader);
		vcn_fem_bcond_destroy(bconditions);
		return NULL;
	}
	/* Allocate data to store Dirichlet Cnd. upon vertices */
	bconditions->Dirichlet_on_vtx_idx = 
		(uint32_t*)malloc(bconditions->N_Dirichlet_on_vtx * sizeof(uint32_t));
	bconditions->Dirichlet_on_vtx_dof_mask =
		(bool*)malloc(2 * bconditions->N_Dirichlet_on_vtx * sizeof(bool));
	bconditions->Dirichlet_on_vtx_val =
		(double*)calloc(2 * bconditions->N_Dirichlet_on_vtx, sizeof(double));

	/* Read number of Neuman conditions upon vertices */
	if(vcn_cfreader_read_uint32_t(cfreader, &(bconditions->N_Neuman_on_vtx)) != 0){
		vcn_cfreader_destroy(cfreader);
		vcn_fem_bcond_destroy(bconditions);
		return NULL;
	}
	/* Allocate data to store Neuman Cnd. upon vertices */
	bconditions->Neuman_on_vtx_idx = 
		(uint32_t*)malloc(bconditions->N_Neuman_on_vtx * sizeof(uint32_t));
	bconditions->Neuman_on_vtx_dof_mask =
		(bool*)malloc(2 * bconditions->N_Neuman_on_vtx * sizeof(bool));
	bconditions->Neuman_on_vtx_val =
		(double*)calloc(2 * bconditions->N_Neuman_on_vtx, sizeof(double));

	/* Read number of Dirichlet conditions upon segments */
	if(vcn_cfreader_read_uint32_t(cfreader, &(bconditions->N_Dirichlet_on_sgm)) != 0){
		vcn_cfreader_destroy(cfreader);
		vcn_fem_bcond_destroy(bconditions);
		return NULL;
	}
	/* Allocate data to store Dirichlet Cnd. upon segments */
	bconditions->Dirichlet_on_sgm_idx = 
		(uint32_t*)malloc(bconditions->N_Dirichlet_on_sgm * sizeof(uint32_t));
	bconditions->Dirichlet_on_sgm_dof_mask =
		(bool*)malloc(2 * bconditions->N_Dirichlet_on_sgm * sizeof(bool));
	bconditions->Dirichlet_on_sgm_val =
		(double*)calloc(2 * bconditions->N_Dirichlet_on_sgm, sizeof(double));

	/* Read number of Neuman conditions upon segments */
	if(vcn_cfreader_read_uint32_t(cfreader, &(bconditions->N_Neuman_on_sgm)) != 0){
		vcn_cfreader_destroy(cfreader);
		vcn_fem_bcond_destroy(bconditions);
		return NULL;
	}
	/* Allocate data to store Neuman Cnd. upon vertices */
	bconditions->Neuman_on_sgm_idx = 
		(uint32_t*)malloc(bconditions->N_Neuman_on_sgm * sizeof(uint32_t));
	bconditions->Neuman_on_sgm_dof_mask =
		(bool*)malloc(2 * bconditions->N_Neuman_on_sgm * sizeof(bool));
	bconditions->Neuman_on_sgm_val =
		(double*)calloc(2 * bconditions->N_Neuman_on_sgm, sizeof(double));

	/* Read Dirichlet Cnd. upon vertices */
	for(uint32_t i = 0; i < bconditions->N_Dirichlet_on_vtx; i++){
		if(vcn_cfreader_read_uint32_t(cfreader, &(bconditions->Dirichlet_on_vtx_idx[i])) != 0){
			vcn_cfreader_destroy(cfreader);
			vcn_fem_bcond_destroy(bconditions);
			return NULL;
		}
		uint32_t aux;
		if(vcn_cfreader_read_uint32_t(cfreader, &aux) != 0){
			vcn_cfreader_destroy(cfreader);
			vcn_fem_bcond_destroy(bconditions);
			return NULL;
		}
		bconditions->Dirichlet_on_vtx_dof_mask[i * 2] = (aux == 1)?true:false;

		if(vcn_cfreader_read_uint32_t(cfreader, &aux) != 0){
			vcn_cfreader_destroy(cfreader);
			vcn_fem_bcond_destroy(bconditions);
			return NULL;
		}
		bconditions->Dirichlet_on_vtx_dof_mask[i*2+1] = (aux == 1)?true:false;

		if(bconditions->Dirichlet_on_vtx_dof_mask[i * 2]){
			if(vcn_cfreader_read_double(cfreader, &(bconditions->Dirichlet_on_vtx_val[i * 2])) != 0){
				vcn_cfreader_destroy(cfreader);
				vcn_fem_bcond_destroy(bconditions);
				return NULL;
			}
		}

		if(bconditions->Dirichlet_on_vtx_dof_mask[i*2+1]){
			if(vcn_cfreader_read_double(cfreader, &(bconditions->Dirichlet_on_vtx_val[i*2+1])) != 0){
				vcn_cfreader_destroy(cfreader);
				vcn_fem_bcond_destroy(bconditions);
				return NULL;
			}
		}
	}

	/* Read Neuman Cnd. upon vertices */
	for(uint32_t i = 0; i < bconditions->N_Neuman_on_vtx; i++){
		if(vcn_cfreader_read_uint32_t(cfreader, &(bconditions->Neuman_on_vtx_idx[i])) != 0){
			vcn_cfreader_destroy(cfreader);
			vcn_fem_bcond_destroy(bconditions);
			return NULL;
		}
		uint32_t aux;
		if(vcn_cfreader_read_uint32_t(cfreader, &aux) != 0){
			vcn_cfreader_destroy(cfreader);
			vcn_fem_bcond_destroy(bconditions);
			return NULL;
		}
		bconditions->Neuman_on_vtx_dof_mask[i * 2] = (aux == 1)?true:false;

		if(vcn_cfreader_read_uint32_t(cfreader, &aux) != 0){
			vcn_cfreader_destroy(cfreader);
			vcn_fem_bcond_destroy(bconditions);
			return NULL;
		}
		bconditions->Neuman_on_vtx_dof_mask[i*2+1] = (aux == 1)?true:false;

		if(bconditions->Neuman_on_vtx_dof_mask[i * 2]){
			if(vcn_cfreader_read_double(cfreader, &(bconditions->Neuman_on_vtx_val[i * 2])) != 0){
				vcn_cfreader_destroy(cfreader);
				vcn_fem_bcond_destroy(bconditions);
				return NULL;
			}
		}

		if(bconditions->Neuman_on_vtx_dof_mask[i*2+1]){
			if(vcn_cfreader_read_double(cfreader, &(bconditions->Neuman_on_vtx_val[i*2+1])) != 0){
				vcn_cfreader_destroy(cfreader);
				vcn_fem_bcond_destroy(bconditions);
				return NULL;
			}
		}
	}

	/* Read Dirichlet Cnd. upon segments */
	for(uint32_t i = 0; i < bconditions->N_Dirichlet_on_sgm; i++){
		if(vcn_cfreader_read_uint32_t(cfreader, &(bconditions->Dirichlet_on_sgm_idx[i])) != 0){
			vcn_cfreader_destroy(cfreader);
			vcn_fem_bcond_destroy(bconditions);
			return NULL;
		}
		uint32_t aux;
		if(vcn_cfreader_read_uint32_t(cfreader, &aux) != 0){
			vcn_cfreader_destroy(cfreader);
			vcn_fem_bcond_destroy(bconditions);
			return NULL;
		}
		bconditions->Dirichlet_on_sgm_dof_mask[i * 2] = (aux == 1)?true:false;

		if(vcn_cfreader_read_uint32_t(cfreader, &aux) != 0){
			vcn_cfreader_destroy(cfreader);
			vcn_fem_bcond_destroy(bconditions);
			return NULL;
		}
		bconditions->Dirichlet_on_sgm_dof_mask[i*2+1] = (aux == 1)?true:false;

		if(bconditions->Dirichlet_on_sgm_dof_mask[i * 2]){
			if(vcn_cfreader_read_double(cfreader, &(bconditions->Dirichlet_on_sgm_val[i * 2])) != 0){
				vcn_cfreader_destroy(cfreader);
				vcn_fem_bcond_destroy(bconditions);
				return NULL;
			}
		}

		if(bconditions->Dirichlet_on_sgm_dof_mask[i*2+1]){
			if(vcn_cfreader_read_double(cfreader, &(bconditions->Dirichlet_on_sgm_val[i*2+1])) != 0){
				vcn_cfreader_destroy(cfreader);
				vcn_fem_bcond_destroy(bconditions);
				return NULL;
			}
		}
	}

	/* Read Neuman Cnd. upon vertices */
	for(uint32_t i = 0; i < bconditions->N_Neuman_on_sgm; i++){
		if(vcn_cfreader_read_uint32_t(cfreader, &(bconditions->Neuman_on_sgm_idx[i])) != 0){
			vcn_cfreader_destroy(cfreader);
			vcn_fem_bcond_destroy(bconditions);
			return NULL;
		}
		uint32_t aux;
		if(vcn_cfreader_read_uint32_t(cfreader, &aux) != 0){
			vcn_cfreader_destroy(cfreader);
			vcn_fem_bcond_destroy(bconditions);
			return NULL;
		}
		bconditions->Neuman_on_sgm_dof_mask[i * 2] = (aux == 1)?true:false;

		if(vcn_cfreader_read_uint32_t(cfreader, &aux) != 0){
			vcn_cfreader_destroy(cfreader);
			vcn_fem_bcond_destroy(bconditions);
			return NULL;
		}
		bconditions->Neuman_on_sgm_dof_mask[i*2+1] = (aux == 1)?true:false;

		if(bconditions->Neuman_on_sgm_dof_mask[i * 2]){
			if(vcn_cfreader_read_double(cfreader, &(bconditions->Neuman_on_sgm_val[i * 2])) != 0){
				vcn_cfreader_destroy(cfreader);
				vcn_fem_bcond_destroy(bconditions);
				return NULL;
			}
		}

		if(bconditions->Neuman_on_sgm_dof_mask[i*2+1]){
			if(vcn_cfreader_read_double(cfreader, &(bconditions->Neuman_on_sgm_val[i*2+1])) != 0){
				vcn_cfreader_destroy(cfreader);
				vcn_fem_bcond_destroy(bconditions);
				return NULL;
			}
		}
	}

	/* Free memory */
	vcn_cfreader_destroy(cfreader);

	/* Return bconditions */
	return bconditions;
}

void vcn_fem_bcond_save(const vcn_bcond_t *const bconditions, 
			const char* filename){
	FILE* fp = fopen(filename, "w");
	if(fp == NULL) return;
  
	fprintf(fp, "# Boundary conditions description \n");
	fprintf(fp, "# Row <- [id x_mask y_mask (x_value) (y_value)]\n\n");

	fprintf(fp, "%i # Number of Dirichlet conditions upon vertices\n",
		bconditions->N_Dirichlet_on_vtx);
	fprintf(fp, "%i # Number of Neuman conditions upon vertices\n",
		bconditions->N_Neuman_on_vtx);
	fprintf(fp, "%i # Number of Dirichlet conditions upon segments\n",
		bconditions->N_Dirichlet_on_sgm);
	fprintf(fp, "%i # Number of Neuman conditions upon segments\n\n",
		bconditions->N_Neuman_on_sgm);

	fprintf(fp, "# Dirichlet conditions upon vertices  \n");
	for(uint32_t i = 0; i < bconditions->N_Dirichlet_on_vtx; i++){
		fprintf(fp, "%i %i %i ",
			bconditions->Dirichlet_on_vtx_idx[i],
			((bconditions->Dirichlet_on_vtx_dof_mask[i * 2])?1:0),
			((bconditions->Dirichlet_on_vtx_dof_mask[i*2+1])?1:0));
		if(bconditions->Dirichlet_on_vtx_dof_mask[i * 2])
			fprintf(fp, "%e ", bconditions->Dirichlet_on_vtx_val[i * 2]);
		if(bconditions->Dirichlet_on_vtx_dof_mask[i*2+1])
			fprintf(fp, "%e ", bconditions->Dirichlet_on_vtx_val[i*2+1]);
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");

	fprintf(fp, "# Neuman conditions upon vertices  \n");
	for(uint32_t i = 0; i < bconditions->N_Neuman_on_vtx; i++){
		fprintf(fp, "%i %i %i ",
			bconditions->Neuman_on_vtx_idx[i],
			((bconditions->Neuman_on_vtx_dof_mask[i * 2])?1:0),
			((bconditions->Neuman_on_vtx_dof_mask[i*2+1])?1:0));
		if(bconditions->Neuman_on_vtx_dof_mask[i * 2])
			fprintf(fp, "%e ", bconditions->Neuman_on_vtx_val[i * 2]);
		if(bconditions->Neuman_on_vtx_dof_mask[i*2+1])
			fprintf(fp, "%e ", bconditions->Neuman_on_vtx_val[i*2+1]);
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");

	fprintf(fp, "# Dirichlet conditions upon segments  \n");
	for(uint32_t i = 0; i < bconditions->N_Dirichlet_on_sgm; i++){
		fprintf(fp, "%i %i %i ",
			bconditions->Dirichlet_on_sgm_idx[i],
			((bconditions->Dirichlet_on_sgm_dof_mask[i * 2])?1:0),
			((bconditions->Dirichlet_on_sgm_dof_mask[i*2+1])?1:0));
		if(bconditions->Dirichlet_on_sgm_dof_mask[i * 2])
			fprintf(fp, "%e ", bconditions->Dirichlet_on_sgm_val[i * 2]);
		if(bconditions->Dirichlet_on_sgm_dof_mask[i*2+1])
			fprintf(fp, "%e ", bconditions->Dirichlet_on_sgm_val[i*2+1]);
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");

	fprintf(fp, "# Neuman conditions upon segments  \n");
	for(uint32_t i = 0; i < bconditions->N_Neuman_on_sgm; i++){
		fprintf(fp, "%i %i %i ",
			bconditions->Neuman_on_sgm_idx[i],
			((bconditions->Neuman_on_sgm_dof_mask[i * 2])?1:0),
			((bconditions->Neuman_on_sgm_dof_mask[i*2+1])?1:0));
		if(bconditions->Neuman_on_sgm_dof_mask[i * 2])
			fprintf(fp, "%e ", bconditions->Neuman_on_sgm_val[i * 2]);
		if(bconditions->Neuman_on_sgm_dof_mask[i*2+1])
			fprintf(fp, "%e ", bconditions->Neuman_on_sgm_val[i*2+1]);
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");  

	fclose(fp);
}

void vcn_fem_bcond_copy(vcn_bcond_t* cpy,
			const vcn_bcond_t *const src){
	cpy->N_dof = src->N_dof;
	cpy->N_Dirichlet_on_vtx = src->N_Dirichlet_on_vtx;
	cpy->N_Neuman_on_vtx = src->N_Neuman_on_vtx;
	cpy->N_Dirichlet_on_sgm = src->N_Dirichlet_on_sgm;
	cpy->N_Neuman_on_sgm = src->N_Neuman_on_sgm;

	if(cpy->N_Dirichlet_on_vtx > 0){
		cpy->Dirichlet_on_vtx_idx =
			(uint32_t*) malloc(cpy->N_Dirichlet_on_vtx * sizeof(uint32_t));
		cpy->Dirichlet_on_vtx_dof_mask =
			(bool*) malloc(cpy->N_dof * cpy->N_Dirichlet_on_vtx * sizeof(bool));
		cpy->Dirichlet_on_vtx_val =
			(double*) malloc(cpy->N_dof * cpy->N_Dirichlet_on_vtx * sizeof(double));

		memcpy(cpy->Dirichlet_on_vtx_idx, src->Dirichlet_on_vtx_idx,
		       cpy->N_Dirichlet_on_vtx * sizeof(uint32_t));
		memcpy(cpy->Dirichlet_on_vtx_dof_mask, src->Dirichlet_on_vtx_dof_mask,
		       cpy->N_dof * cpy->N_Dirichlet_on_vtx * sizeof(bool));
		memcpy(cpy->Dirichlet_on_vtx_val, src->Dirichlet_on_vtx_val,
		       cpy->N_dof * cpy->N_Dirichlet_on_vtx * sizeof(double));
	}

	if(cpy->N_Neuman_on_vtx > 0){
		cpy->Neuman_on_vtx_idx =
			(uint32_t*) malloc(cpy->N_Neuman_on_vtx * sizeof(uint32_t));
		cpy->Neuman_on_vtx_dof_mask =
			(bool*) malloc(cpy->N_dof * cpy->N_Neuman_on_vtx * sizeof(bool));
		cpy->Neuman_on_vtx_val =
			(double*) malloc(cpy->N_dof * cpy->N_Neuman_on_vtx * sizeof(double));

		memcpy(cpy->Neuman_on_vtx_idx, src->Neuman_on_vtx_idx,
		       cpy->N_Neuman_on_vtx * sizeof(uint32_t));
		memcpy(cpy->Neuman_on_vtx_dof_mask, src->Neuman_on_vtx_dof_mask,
		       cpy->N_dof * cpy->N_Neuman_on_vtx * sizeof(bool));
		memcpy(cpy->Neuman_on_vtx_val, src->Neuman_on_vtx_val,
		       cpy->N_dof * cpy->N_Neuman_on_vtx * sizeof(double));
	}

	if(cpy->N_Dirichlet_on_sgm > 0){
		cpy->Dirichlet_on_sgm_idx =
			(uint32_t*) malloc(cpy->N_Dirichlet_on_sgm * sizeof(uint32_t));
		cpy->Dirichlet_on_sgm_dof_mask =
			(bool*) malloc(cpy->N_dof * cpy->N_Dirichlet_on_sgm * sizeof(bool));
		cpy->Dirichlet_on_sgm_val =
			(double*) malloc(cpy->N_dof * cpy->N_Dirichlet_on_sgm * sizeof(double));

		memcpy(cpy->Dirichlet_on_sgm_idx, src->Dirichlet_on_sgm_idx,
		       cpy->N_Dirichlet_on_sgm * sizeof(uint32_t));
		memcpy(cpy->Dirichlet_on_sgm_dof_mask, src->Dirichlet_on_sgm_dof_mask,
		       cpy->N_dof * cpy->N_Dirichlet_on_sgm * sizeof(bool));
		memcpy(cpy->Dirichlet_on_sgm_val, src->Dirichlet_on_sgm_val,
		       cpy->N_dof * cpy->N_Dirichlet_on_sgm * sizeof(double));
	}

	if(cpy->N_Neuman_on_sgm > 0){
		cpy->Neuman_on_sgm_idx =
			(uint32_t*) malloc(cpy->N_Neuman_on_sgm * sizeof(uint32_t));
		cpy->Neuman_on_sgm_dof_mask =
			(bool*) malloc(cpy->N_dof * cpy->N_Neuman_on_sgm * sizeof(bool));
		cpy->Neuman_on_sgm_val =
			(double*) malloc(cpy->N_dof * cpy->N_Neuman_on_sgm * sizeof(double));

		memcpy(cpy->Neuman_on_sgm_idx, src->Neuman_on_sgm_idx,
		       cpy->N_Neuman_on_sgm * sizeof(uint32_t));
		memcpy(cpy->Neuman_on_sgm_dof_mask, src->Neuman_on_sgm_dof_mask,
		       cpy->N_dof * cpy->N_Neuman_on_sgm * sizeof(bool));
		memcpy(cpy->Neuman_on_sgm_val, src->Neuman_on_sgm_val,
		       cpy->N_dof * cpy->N_Neuman_on_sgm * sizeof(double));
	}
}

void vcn_fem_bcond_clear(vcn_bcond_t* bconditions){
	if(bconditions->N_Dirichlet_on_vtx > 0){
		free(bconditions->Dirichlet_on_vtx_idx);
		free(bconditions->Dirichlet_on_vtx_dof_mask);
		free(bconditions->Dirichlet_on_vtx_val);
	}
	if(bconditions->N_Neuman_on_vtx > 0){
		free(bconditions->Neuman_on_vtx_idx);
		free(bconditions->Neuman_on_vtx_dof_mask);
		free(bconditions->Neuman_on_vtx_val);
	}
	if(bconditions->N_Dirichlet_on_sgm > 0){
		free(bconditions->Dirichlet_on_sgm_idx);
		free(bconditions->Dirichlet_on_sgm_dof_mask);
		free(bconditions->Dirichlet_on_sgm_val);
	}
	if(bconditions->N_Neuman_on_sgm > 0){
		free(bconditions->Neuman_on_sgm_idx);
		free(bconditions->Neuman_on_sgm_dof_mask);
		free(bconditions->Neuman_on_sgm_val);
	}
	memset(bconditions, 0, sizeof(vcn_bcond_t));
}

void vcn_fem_bcond_destroy(vcn_bcond_t* bconditions){
	vcn_fem_bcond_clear(bconditions);
	free(bconditions);
}

void vcn_fem_bcond_printf(const vcn_bcond_t* const bconditions){
	printf("Dirichlet conditions on vertices: %i \n", bconditions->N_Dirichlet_on_vtx);
	for(uint32_t i = 0; i < bconditions->N_Dirichlet_on_vtx; i++){
		printf("  %i (%i %i) <- %e %e\n",
		       bconditions->Dirichlet_on_vtx_idx[i],
		       bconditions->Dirichlet_on_vtx_dof_mask[i * 2],
		       bconditions->Dirichlet_on_vtx_dof_mask[i*2+1],
		       bconditions->Dirichlet_on_vtx_val[i * 2],
		       bconditions->Dirichlet_on_vtx_val[i*2+1]);
	}
	printf("Dirichlet conditions on segments: %i \n", bconditions->N_Dirichlet_on_sgm);
	for(uint32_t i = 0; i < bconditions->N_Dirichlet_on_sgm; i++){
		printf("  %i (%i %i) <- %e %e\n",
		       bconditions->Dirichlet_on_sgm_idx[i],
		       bconditions->Dirichlet_on_sgm_dof_mask[i * 2],
		       bconditions->Dirichlet_on_sgm_dof_mask[i*2+1],
		       bconditions->Dirichlet_on_sgm_val[i * 2],
		       bconditions->Dirichlet_on_sgm_val[i*2+1]);
	}
	printf("Neumann conditions on vertices: %i \n", bconditions->N_Neuman_on_vtx);
	for(uint32_t i = 0; i < bconditions->N_Neuman_on_vtx; i++){
		printf("  %i (%i %i) <- %e %e\n",
		       bconditions->Neuman_on_vtx_idx[i],
		       bconditions->Neuman_on_vtx_dof_mask[i * 2],
		       bconditions->Neuman_on_vtx_dof_mask[i*2+1],
		       bconditions->Neuman_on_vtx_val[i * 2],
		       bconditions->Neuman_on_vtx_val[i*2+1]);
	}
	printf("Neuman conditions on segments: %i \n", bconditions->N_Neuman_on_sgm);
	for(uint32_t i = 0; i < bconditions->N_Neuman_on_sgm; i++){
		printf("  %i (%i %i) <- %e %e\n",
		       bconditions->Neuman_on_sgm_idx[i],
		       bconditions->Neuman_on_sgm_dof_mask[i * 2],
		       bconditions->Neuman_on_sgm_dof_mask[i*2+1],
		       bconditions->Neuman_on_sgm_val[i * 2],
		       bconditions->Neuman_on_sgm_val[i*2+1]);
	}
}

vcn_bcond_t* vcn_fem_bcond_create_from_model_to_mesh
(const vcn_msh3trg_t *const msh3trg,
 const vcn_bcond_t *const bcond_on_model)
/* Mapping of the boundary conditions on the model to the boundary 
 * conditions on the mesh. The conditions are placed on the nodes.
 */
{
	vcn_bcond_t* bmeshcond = 
		(vcn_bcond_t*) calloc(1, sizeof(vcn_bcond_t));

	bmeshcond->N_dof = bcond_on_model->N_dof;
	/***** Get Dirichlet conditions *****/
	bmeshcond->N_Dirichlet_on_vtx = bcond_on_model->N_Dirichlet_on_vtx;
	/* Set Dirichlet conditions on vertices forming the segments with 
	 * the original conditions.
	 */
	for(uint32_t i=0; i < bcond_on_model->N_Dirichlet_on_sgm; i++){
		uint32_t geom_idx = bcond_on_model->Dirichlet_on_sgm_idx[i];
		if(msh3trg->N_subsgm_x_inputsgm[geom_idx] == 0) continue;
		bmeshcond->N_Dirichlet_on_vtx += 
			msh3trg->N_subsgm_x_inputsgm[geom_idx] + 1;
	}
	/* Allocate Dirichlet conditions */
	if(bmeshcond->N_Dirichlet_on_vtx > 0){
		bmeshcond->Dirichlet_on_vtx_idx =
			(uint32_t*)calloc(bmeshcond->N_Dirichlet_on_vtx,
					  sizeof(uint32_t));
		bmeshcond->Dirichlet_on_vtx_dof_mask =
			(bool*)malloc(bmeshcond->N_dof * bmeshcond->N_Dirichlet_on_vtx *
				      sizeof(bool));
		bmeshcond->Dirichlet_on_vtx_val =
			(double*)malloc(bmeshcond->N_dof * bmeshcond->N_Dirichlet_on_vtx *
					sizeof(double));

		/* Set Dirichlet conditions on vertices */
		for(uint32_t i = 0; i < bcond_on_model->N_Dirichlet_on_vtx; i++){
			uint32_t geom_idx = bcond_on_model->Dirichlet_on_vtx_idx[i];
			bmeshcond->Dirichlet_on_vtx_idx[i] =
				msh3trg->input_vertices[geom_idx];
		}

		memcpy(bmeshcond->Dirichlet_on_vtx_dof_mask,
		       bcond_on_model->Dirichlet_on_vtx_dof_mask, sizeof(bool) *
		       bcond_on_model->N_Dirichlet_on_vtx * bcond_on_model->N_dof);
		memcpy(bmeshcond->Dirichlet_on_vtx_val,
		       bcond_on_model->Dirichlet_on_vtx_val, sizeof(double) *
		       bcond_on_model->N_Dirichlet_on_vtx * bcond_on_model->N_dof);

		/* Set Dirichlet conditions on segments */
		uint32_t idx = bcond_on_model->N_Dirichlet_on_vtx;
		for(uint32_t i=0; i < bcond_on_model->N_Dirichlet_on_sgm; i++){
			uint32_t geom_idx = bcond_on_model->Dirichlet_on_sgm_idx[i];
			if(msh3trg->N_subsgm_x_inputsgm[geom_idx] == 0) continue;
			for(uint32_t j = 0; j < msh3trg->N_subsgm_x_inputsgm[geom_idx] + 1; j++){
				bmeshcond->Dirichlet_on_vtx_idx[idx] = 
					msh3trg->meshvtx_x_inputsgm[geom_idx][j];

				for(uint32_t k = 0; k < bcond_on_model->N_dof; k++){
					bmeshcond->Dirichlet_on_vtx_dof_mask[idx * bcond_on_model->N_dof + k] = 
						bcond_on_model->Dirichlet_on_sgm_dof_mask[i * bcond_on_model->N_dof + k];
				}

				for(uint32_t k = 0; k < bcond_on_model->N_dof; k++){	
					bmeshcond->Dirichlet_on_vtx_val[idx * bcond_on_model->N_dof + k] = 
						bcond_on_model->Dirichlet_on_sgm_val[i * bcond_on_model->N_dof + k];
				}

				idx ++;
			}
		}
	}

	/***** Get Neuman conditions *****/
	bmeshcond->N_Neuman_on_vtx = bcond_on_model->N_Neuman_on_vtx;
	/* Set Neuman conditions on vertices forming the segments with 
	 * the original conditions.
	 */
	for(uint32_t i=0; i < bcond_on_model->N_Neuman_on_sgm; i++){
		uint32_t geom_idx = bcond_on_model->Neuman_on_sgm_idx[i];
		if(msh3trg->N_subsgm_x_inputsgm[geom_idx] == 0) continue;
		bmeshcond->N_Neuman_on_vtx += 
			msh3trg->N_subsgm_x_inputsgm[geom_idx] + 1;
	}
	/* Allocate Neuman conditions */
	if(bmeshcond->N_Neuman_on_vtx > 0){
		bmeshcond->Neuman_on_vtx_idx =
			(uint32_t*)calloc(bmeshcond->N_Neuman_on_vtx, sizeof(uint32_t));
		bmeshcond->Neuman_on_vtx_dof_mask =
			(bool*)malloc(bmeshcond->N_dof * bmeshcond->N_Neuman_on_vtx *
				      sizeof(bool));
		bmeshcond->Neuman_on_vtx_val =
			(double*)malloc(bmeshcond->N_dof * bmeshcond->N_Neuman_on_vtx *
					sizeof(double));
    
		/* Set Neuman conditions on vertices */
		for(uint32_t i=0; i < bcond_on_model->N_Neuman_on_vtx; i++){
			uint32_t geom_idx = bcond_on_model->Neuman_on_vtx_idx[i];
			bmeshcond->Neuman_on_vtx_idx[i] =
				msh3trg->input_vertices[geom_idx];
		}

		memcpy(bmeshcond->Neuman_on_vtx_dof_mask,
		       bcond_on_model->Neuman_on_vtx_dof_mask, sizeof(bool) *
		       bcond_on_model->N_Neuman_on_vtx * bcond_on_model->N_dof);
		memcpy(bmeshcond->Neuman_on_vtx_val,
		       bcond_on_model->Neuman_on_vtx_val, sizeof(double) *
		       bcond_on_model->N_Neuman_on_vtx * bcond_on_model->N_dof);

		/* Set Neuman conditions on segments */
		uint32_t idx = bcond_on_model->N_Neuman_on_vtx;
		for(uint32_t i=0; i < bcond_on_model->N_Neuman_on_sgm; i++){
			uint32_t geom_idx = bcond_on_model->Neuman_on_sgm_idx[i];
			if(msh3trg->N_subsgm_x_inputsgm[geom_idx] == 0) continue;
      
			uint32_t s1 = msh3trg->meshvtx_x_inputsgm[geom_idx][0];
			uint32_t s2 = msh3trg->meshvtx_x_inputsgm[geom_idx]
				[msh3trg->N_subsgm_x_inputsgm[geom_idx]];
			double sgm_length = vcn_utils2D_get_dist(&(msh3trg->vertices[s1 * 2]),
								 &(msh3trg->vertices[s2 * 2]));

			for(uint32_t j = 0; j < msh3trg->N_subsgm_x_inputsgm[geom_idx] + 1; j++){
				uint32_t sj = msh3trg->meshvtx_x_inputsgm[geom_idx][j];
				bmeshcond->Neuman_on_vtx_idx[idx] = sj;

				for(uint32_t k = 0; k < bcond_on_model->N_dof; k++){
					bmeshcond->Neuman_on_vtx_dof_mask[idx * bcond_on_model->N_dof + k] = 
						bcond_on_model->Neuman_on_sgm_dof_mask[i * bcond_on_model->N_dof + k];
				}

				/* Compute proportional weight to assign condition */
				double w_left = 0.0;
				if(j > 0){
					uint32_t sj_prev = msh3trg->meshvtx_x_inputsgm[geom_idx][j-1];
					double dist =
						vcn_utils2D_get_dist(&(msh3trg->vertices[sj_prev * 2]),
								     &(msh3trg->vertices[sj * 2]));
					w_left = 0.5 * (dist);
				}

				double w_right = 0.0;
				if(j < msh3trg->N_subsgm_x_inputsgm[geom_idx]){
					uint32_t sj_next = msh3trg->meshvtx_x_inputsgm[geom_idx][j+1];
					double dist =
						vcn_utils2D_get_dist(&(msh3trg->vertices[sj_next * 2]),
								     &(msh3trg->vertices[sj * 2]));
					w_right = 0.5 * (dist);
				}
	
				double weight = (w_left + w_right) / sgm_length;

				for(uint32_t k = 0; k < bcond_on_model->N_dof; k++){	
					bmeshcond->Neuman_on_vtx_val[idx * bcond_on_model->N_dof + k] = 
						bcond_on_model->Neuman_on_sgm_val[i * bcond_on_model->N_dof + k] * weight;
				}

				idx ++;
			}
		}
	}

	return bmeshcond;
}
  
vcn_fem_elem_t* vcn_fem_elem_create_triangle(){
	/* Triangular element of linear interpolation */
	vcn_fem_elem_t* elemtype = (vcn_fem_elem_t*)
		malloc(sizeof(vcn_fem_elem_t));
	elemtype->N_nodes = 3;
	elemtype->N_Gauss_points = 1;

	elemtype->Ni = 
		(double (**)(double, double))calloc(elemtype->N_nodes, sizeof(void*));
	elemtype->dNi_dpsi = 
		(double (**)(double, double))calloc(elemtype->N_nodes, sizeof(void*));
	elemtype->dNi_deta = 
		(double (**)(double, double))calloc(elemtype->N_nodes, sizeof(void*));  
	elemtype->Ni[0] = N1;
	elemtype->Ni[1] = N2;
	elemtype->Ni[2] = N3;
	elemtype->dNi_dpsi[0] = dN1_dpsi;
	elemtype->dNi_dpsi[1] = dN2_dpsi;
	elemtype->dNi_dpsi[2] = dN3_dpsi;
	elemtype->dNi_deta[0] = dN1_deta;
	elemtype->dNi_deta[1] = dN2_deta;
	elemtype->dNi_deta[2] = dN3_deta;
  
	elemtype->psi = 
		(double*)malloc(elemtype->N_Gauss_points*sizeof(double));
	elemtype->eta = 
		(double*)malloc(elemtype->N_Gauss_points*sizeof(double));
	elemtype->gp_weight = 
		(double*)malloc(elemtype->N_Gauss_points*sizeof(double));

	elemtype->psi[0] = 1.0/3.0;
	elemtype->eta[0] = 1.0/3.0;

	elemtype->gp_weight[0] = 0.5;

	return elemtype;
}

void vcn_fem_elem_destroy(vcn_fem_elem_t* elemtype){
	free(elemtype->Ni);
	free(elemtype->dNi_dpsi);
	free(elemtype->dNi_deta);
	free(elemtype->psi);
	free(elemtype->eta);
	free(elemtype->gp_weight);
	free(elemtype);
}

uint32_t vcn_fem_elem_get_N_nodes(const vcn_fem_elem_t *const elemtype){
	return elemtype->N_nodes;
}

void* vcn_fem_elem_get_ith_shape_function
(const vcn_fem_elem_t *const elemtype, uint32_t i){
	return elemtype->Ni[i];
}

uint32_t vcn_fem_elem_get_closest_Gauss_Point_to_the_ith_node
(const vcn_fem_elem_t *const elemtype, uint32_t i){
	if(elemtype->N_nodes == 3 && 
	   elemtype->N_Gauss_points == 1){
		/* Triangle of linear interpolation
		 *              o
		 *             / \
		 *            / * \    <---- Gauss Point
		 *           /_____\
		 *          o       o  <---- Nodes
		 */
		return 0;
	}
	printf("FEM Error: Unknown element type with %i nodes and %i Gauss points.\n",
	       elemtype->N_nodes, elemtype->N_Gauss_points);
	return 0;
}

static void pipeline_assemble_system
(vcn_sparse_t* K, double* M, double *F,
 const vcn_msh3trg_t *const mesh,
 const vcn_fem_elem_t *const elemtype,
 const vcn_fem_material_t *const material,
 bool enable_self_weight,
 double gravity[2],
 bool enable_plane_stress,
 double thickness,
 bool enable_computing_damage,
 double* damage_elem,
 bool* elements_enabled /* NULL to enable all */){
	uint32_t N_elements = mesh->N_triangles;

	vcn_sparse_reset(K);
	if(M !=NULL) memset(M, 0, vcn_sparse_get_size(K) * sizeof(double));
	memset(F, 0, vcn_sparse_get_size(K) * sizeof(double));
	/* Allocate elemental Stiffness Matrix and Force Vector */
	double* Ke = (double*)malloc(4 * elemtype->N_nodes * elemtype->N_nodes *
				     sizeof(double));
	double* Me = NULL;
	if(M != NULL)
		Me = (double*)malloc(2 * elemtype->N_nodes*sizeof(double));
	double* Fe = (double*)
		malloc(2 * elemtype->N_nodes*sizeof(double));

	if (!enable_plane_stress)
		thickness = 1.0;

	/* Assembly global system */
	uint32_t N_negative_jacobians = 0;
	for(uint32_t k=0; k < N_elements; k++){
		/* Get material properties */
		double density = vcn_fem_material_get_density(material);
		double E = vcn_fem_material_get_elasticity_module(material);

		if(elements_enabled != NULL){
			if(!elements_enabled[k]){
				E *= 1e-6;
				density *= 1e-6;
			}
		}

		double v = vcn_fem_material_get_poisson_module(material);
    
		/* Get constitutive matrix */
		double d11 = E/(1.0 - vcn_math_pow2(v));
		double d12 = v*d11;
    
		if(!enable_plane_stress){
			d11 = (E*(1.0-v))/((1.0 + v)*(1.0-2*v));
			d12 = (v*d11)/(1.0-v);
		}    
		double d22 = d11;
		double d33 = E/(2.0*(1.0+v));

		/* Allocate Cartesian derivatives for each Gauss Point */
		double* dNi_dx = (double*)malloc(elemtype->N_nodes*sizeof(double));
		double* dNi_dy = (double*)malloc(elemtype->N_nodes*sizeof(double));
		/* Compute constitutive matrix */
		double fx = 0.0;
		double fy = 0.0;
		if(enable_self_weight){
			fx = gravity[0] * vcn_fem_material_get_density(material);
			fy = gravity[1] * vcn_fem_material_get_density(material);
		}
    
		/* Integrate Ke and Fe using Gauss quadrature */
		memset(Ke, 0, 4 * elemtype->N_nodes * elemtype->N_nodes *
		       sizeof(double));
		if(M != NULL) 
			memset(Me, 0, 2 * elemtype->N_nodes * sizeof(double));
		memset(Fe, 0, 2 * elemtype->N_nodes * sizeof(double));
		for(uint32_t j=0; j < elemtype->N_Gauss_points; j++){      
			/* Get constitutive model */
			double dd11  = d11;
			double dd12  = d12;
			double dd22  = d22;
			double dd33  = d33;
      
			if(enable_computing_damage){
				dd11 *= (1.0 - damage_elem[k*elemtype->N_Gauss_points + j]);
				dd12 *= (1.0 - damage_elem[k*elemtype->N_Gauss_points + j]);
				dd22 *= (1.0 - damage_elem[k*elemtype->N_Gauss_points + j]);
				dd33 *= (1.0 - damage_elem[k*elemtype->N_Gauss_points + j]);
			}

			/* Compute Jacobian derivatives */
			double dx_dpsi = 0.0;
			double dy_dpsi = 0.0;
			double dx_deta = 0.0;
			double dy_deta = 0.0;
			for(uint32_t i=0; i < elemtype->N_nodes; i++){
				uint32_t inode = mesh->vertices_forming_triangles[k * elemtype->N_nodes + i];
				double xi = mesh->vertices[inode * 2];
				double yi = mesh->vertices[inode*2+1];
				dx_dpsi += 
					elemtype->dNi_dpsi[i](elemtype->psi[j], elemtype->eta[j]) * xi;
				dx_deta +=
					elemtype->dNi_deta[i](elemtype->psi[j], elemtype->eta[j]) * xi;
				dy_dpsi += 
					elemtype->dNi_dpsi[i](elemtype->psi[j], elemtype->eta[j]) * yi;
				dy_deta +=
					elemtype->dNi_deta[i](elemtype->psi[j], elemtype->eta[j]) * yi;
			}

			/* Compute Jacobian inverse and determinant */
			double detJ = dx_dpsi*dy_deta - dy_dpsi*dx_deta;
			/* Check if the element is distorted */
			if(detJ < 0) N_negative_jacobians += 1;
      
			double Jinv[4];
			Jinv[0] =  dy_deta/detJ;
			Jinv[1] = -dy_dpsi/detJ;
			Jinv[2] = -dx_deta/detJ;
			Jinv[3] =  dx_dpsi/detJ;
      
			/* Compute Shape functions derivatives in cartesian space */ 
			for(uint32_t i=0; i < elemtype->N_nodes; i++){
				dNi_dx[i] = 
					Jinv[0]*elemtype->dNi_dpsi[i](elemtype->psi[j], elemtype->eta[j]) + 
					Jinv[1]*elemtype->dNi_deta[i](elemtype->psi[j], elemtype->eta[j]);
				dNi_dy[i] = 
					Jinv[2]*elemtype->dNi_dpsi[i](elemtype->psi[j], elemtype->eta[j]) + 
					Jinv[3]*elemtype->dNi_deta[i](elemtype->psi[j], elemtype->eta[j]);
			}
			/* Compute elemental stiffness matrix
			 *       _            _         _           _        _  _
			 *      |dNi/dx        |       | d11 d12     |      |  |
			 * Bi = |       dNi/dy |   D = | d21 d22     |  K = |  | B'DB dx dy
			 *      |dNi/dy dNi/dx_|       |_        d33_|     _| _|

			 * B  = [B1 B2...Bi... Belemtype->N_nodes]
			 *
			 */
      
			for(uint32_t i1 = 0; i1 < elemtype->N_nodes; i1++){	
				for(uint32_t i2 = 0; i2 < elemtype->N_nodes; i2++){
					/*  Integrating elemental siffness matrix */
					Ke[(i1 * 2)*(2*elemtype->N_nodes) + (i2 * 2)] += 
						(dNi_dx[i1]*dNi_dx[i2]*dd11 + dNi_dy[i1]*dNi_dy[i2]*dd33) *
						detJ * thickness * elemtype->gp_weight[j];
					Ke[(i1 * 2)*(2*elemtype->N_nodes) + (i2*2+1)] +=
						(dNi_dx[i1]*dNi_dy[i2]*dd12 + dNi_dy[i1]*dNi_dx[i2]*dd33) *
						detJ * thickness * elemtype->gp_weight[j];
					Ke[(i1*2+1)*(2*elemtype->N_nodes) + (i2 * 2)] +=
						(dNi_dy[i1]*dNi_dx[i2]*dd12 + dNi_dx[i1]*dNi_dy[i2]*dd33) *
						detJ * thickness * elemtype->gp_weight[j];
					Ke[(i1*2+1)*(2*elemtype->N_nodes) + (i2*2+1)] +=
						(dNi_dy[i1]*dNi_dy[i2]*dd22 + dNi_dx[i1]*dNi_dx[i2]*dd33) *
						detJ * thickness * elemtype->gp_weight[j];
				}
				/* Calculate shape function of the i-th node at the j-th gauss point */
				double Ni_eval = 
					elemtype->Ni[i1](elemtype->psi[j], elemtype->eta[j]);
				/*  Integrating elemental mass matrix */
				if(M != NULL){
					/* OPPORTUNITY: Allocate just one for each component */
					Me[i1 * 2] += Ni_eval * Ni_eval * density *
						detJ * thickness * elemtype->gp_weight[j];
					Me[i1*2+1] += Ni_eval * Ni_eval * density * 
						detJ * thickness * elemtype->gp_weight[j];
				}
				/* Compute elemental forces vector */
				/* OPPORTUNITY: Allocate just one for each component */
				Fe[i1 * 2] += Ni_eval * fx * detJ * 
					thickness * elemtype->gp_weight[j];
				Fe[i1*2+1] += Ni_eval * fy * detJ * 
					thickness * elemtype->gp_weight[j];
			}
		}
		/* Add to global stiffness matrix */
		for(uint32_t i1 = 0; i1 < elemtype->N_nodes; i1++){
			for(uint32_t i2 = 0; i2 < elemtype->N_nodes; i2++){

				for(uint8_t j1=0; j1 < 2; j1++){
					for(uint8_t j2=0; j2 < 2; j2++){
						vcn_sparse_add
							(K, mesh->vertices_forming_triangles[k*3+i1]*2 + j1,
							 mesh->vertices_forming_triangles[k*3+i2]*2 + j2,
							 Ke[(i1*2+j1)*(2*elemtype->N_nodes) + (i2*2+j2)]);
					}
				}
			}
			/* Add to global mass matrix */
			if (M != NULL) {
				M[mesh->vertices_forming_triangles[k*3+i1] * 2] += Me[i1 * 2];
				M[mesh->vertices_forming_triangles[k*3+i1]*2+1] += Me[i1*2+1];
			}
			/* Add to global forces vector */
			F[mesh->vertices_forming_triangles[k*3+i1] * 2] += Fe[i1 * 2];
			F[mesh->vertices_forming_triangles[k*3+i1]*2+1] += Fe[i1*2+1];
		}
		free(dNi_dx);
		free(dNi_dy);
	}
	if (N_negative_jacobians > 0)
		/* OPPORTUNITY: Send to logfile */
		printf("ERROR: There are %i negative Jacobians.\n", N_negative_jacobians);

	/* Free elemental stiffness matrix and force vector */
	free(Ke);
	if (M != NULL) 
		free(Me);
	free(Fe);
}

static void pipeline_set_boundary_conditions
(vcn_sparse_t* K,
 double* F, 
 const vcn_bcond_t *const bmeshcond, 
 double thickness,
 double factor)
/* We assume that all the conditions are placed on the mesh nodes */
{
	/* OPPORTUNITY: Use libre_evaluator to eval time functions */
	/* Set Neuman conditions on vertices (Loads in the Solid Mechanics context) */
	for(uint32_t i=0; i < bmeshcond->N_Neuman_on_vtx; i++){
		for(uint32_t j=0; j < bmeshcond->N_dof; j++){
			if(!bmeshcond->Neuman_on_vtx_dof_mask[i * bmeshcond->N_dof + j])
				continue;
			F[bmeshcond->Neuman_on_vtx_idx[i] * bmeshcond->N_dof + j] += factor *
				bmeshcond->Neuman_on_vtx_val[i * bmeshcond->N_dof + j];
		}
	}

	/* Set Dirichlet conditions (Displacements in the Solid Mechanics context)  */
	for(uint32_t i=0; i < bmeshcond->N_Dirichlet_on_vtx; i++){
		for(uint32_t k=0; k < bmeshcond->N_dof; k++){
			if(!bmeshcond->Dirichlet_on_vtx_dof_mask[i * bmeshcond->N_dof + k])
				continue;
			uint32_t idx = bmeshcond->Dirichlet_on_vtx_idx[i] * bmeshcond->N_dof + k;
			double value = factor *
				bmeshcond->Dirichlet_on_vtx_val[i * bmeshcond->N_dof + k];
			vcn_sparse_set_Dirichlet_condition(K, F, idx, value);
		}
	}
}

static void pipeline_compute_strain
(double *strain,
 const vcn_msh3trg_t *const mesh,
 double *displacement,
 const vcn_fem_elem_t *const elemtype,
 bool enable_computing_damage,
 bool enable_plane_stress,
 const vcn_fem_material_t *const material,
 double *damage,
 double *r_dmg_prev,
 double *r_dmg)
{
	double *vertices = (double*) mesh->vertices;
	uint32_t N_elements = mesh->N_triangles;
	uint32_t *connectivity_mtx = (uint32_t*) mesh->vertices_forming_triangles;

	/* Initialize strains */
	memset(strain, 0, 3 * elemtype->N_Gauss_points * N_elements * sizeof(double));

	/* Iterate over elements to compute strain, stress and damage at nodes */
	for(uint32_t k=0; k < N_elements; k++){
		double* dNi_dx = (double*)malloc(elemtype->N_nodes*sizeof(double));
		double* dNi_dy = (double*)malloc(elemtype->N_nodes*sizeof(double));

		/* Integrate domain */
		for(uint32_t j=0; j < elemtype->N_Gauss_points; j++){
			uint32_t idx = k * elemtype->N_Gauss_points + j;

			/* Compute Jacobian derivatives */
			double dx_dpsi = 0.0;
			double dy_dpsi = 0.0;
			double dx_deta = 0.0;
			double dy_deta = 0.0;

			for(uint32_t i=0; i < elemtype->N_nodes; i++){
				uint32_t inode = connectivity_mtx[k * elemtype->N_nodes + i];
				double xi = vertices[inode * 2];
				double yi = vertices[inode*2+1];
				dx_dpsi += 
					elemtype->dNi_dpsi[i](elemtype->psi[j], elemtype->eta[j]) * xi;
				dx_deta += 
					elemtype->dNi_deta[i](elemtype->psi[j], elemtype->eta[j]) * xi;
				dy_dpsi +=
					elemtype->dNi_dpsi[i](elemtype->psi[j], elemtype->eta[j]) * yi;
				dy_deta +=
					elemtype->dNi_deta[i](elemtype->psi[j], elemtype->eta[j]) * yi;
			}
      
			/* Compute Jacobian inverse and determinant */
			double detJ = dx_dpsi*dy_deta - dy_dpsi*dx_deta;
			double Jinv[4];
			Jinv[0] =  dy_deta/detJ;
			Jinv[1] = -dy_dpsi/detJ;
			Jinv[2] = -dx_deta/detJ;
			Jinv[3] =  dx_dpsi/detJ;
      
			/* Compute Shape functions derivatives in cartesian space */ 
			for(uint32_t i=0; i < elemtype->N_nodes; i++){
				dNi_dx[i] = 
					Jinv[0]*elemtype->dNi_dpsi[i](elemtype->psi[j], elemtype->eta[j]) + 
					Jinv[1]*elemtype->dNi_deta[i](elemtype->psi[j], elemtype->eta[j]);
				dNi_dy[i] = 
					Jinv[2]*elemtype->dNi_dpsi[i](elemtype->psi[j], elemtype->eta[j]) + 
					Jinv[3]*elemtype->dNi_deta[i](elemtype->psi[j], elemtype->eta[j]);
			}

			/* Compute Strain at Gauss Point */
			for(uint32_t i=0; i < elemtype->N_nodes; i++){
				uint32_t inode = connectivity_mtx[k * elemtype->N_nodes + i];
				strain[idx * 3] += dNi_dx[i] * displacement[inode * 2];
				strain[idx*3+1] += dNi_dy[i] * displacement[inode*2+1];
				strain[idx*3+2] += (dNi_dy[i] * displacement[inode * 2] +
						    dNi_dx[i] * displacement[inode*2+1]);
			}
      
			/* Compute damage */
			if(enable_computing_damage){
				uint32_t n1 = connectivity_mtx[k * 3];
				uint32_t n2 = connectivity_mtx[k*3+1];
				uint32_t n3 = connectivity_mtx[k*3+2];
				double characteristic_length_of_fractured_domain =
					sqrt(((vertices[n2 * 2]-vertices[n1 * 2]) *
					      (vertices[n3*2+1]-vertices[n1*2+1]) -
					      (vertices[n2*2+1]-vertices[n1*2+1]) *
					      (vertices[n3 * 2]-vertices[n1 * 2])));
      
				damage[idx] = 
					material->damage(material, 
							 &(strain[idx*3]),
							 &(r_dmg_prev[idx]),
							 &(r_dmg[idx]),
							 characteristic_length_of_fractured_domain,
							 enable_plane_stress);
			}
		}
		free(dNi_dx);
		free(dNi_dy);
	}
}


static void pipeline_compute_main_stress(double *stress, 
					 double *main_stress,
					 uint32_t N_elements,
					 const vcn_fem_elem_t *const elemtype){
	for(uint32_t i=0; i < N_elements; i++){
		for(uint32_t j=0; j < elemtype->N_Gauss_points; j++){
			uint32_t idx = i*elemtype->N_Gauss_points + j;
			double sigma_avg = 0.5 * (stress[idx * 3] + stress[idx*3+1]);
			double R = sqrt(0.25 * 
					vcn_math_pow2(stress[idx * 3] - stress[idx*3+1]) +
					vcn_math_pow2(stress[idx*3+2]));
			main_stress[idx * 2] = sigma_avg + R;
			main_stress[idx*2+1] = sigma_avg - R;
		}
	}
}

static void pipeline_compute_error_on_elements
(double* error,
 const vcn_msh3trg_t *const mesh,
 double *displacement,
 double* strain,
 const vcn_fem_elem_t *const elemtype)
/* Compute elemental error based on strain field */
{
	uint32_t N_vertices = mesh->N_vertices;
	uint32_t N_elements = mesh->N_triangles;


	uint32_t N_gp = elemtype->N_Gauss_points;
	memset(error, 0, N_elements * N_gp * sizeof(double));

	/* Interpolate strain to nodes */
	double* strain_on_nodes = (double*)calloc(3 * N_vertices, sizeof(double));
	vcn_fem_interpolate_from_Gauss_points_to_nodes(mesh, elemtype,
						       3, strain, strain_on_nodes);

	uint32_t *connectivity_mtx = (uint32_t*) mesh->vertices_forming_triangles;

	/* Return values to Gauss Points to compute the error in the element */
	for(uint32_t i=0; i < N_elements; i++){
		error[i] = 0;
		for(uint32_t j=0; j < N_gp; j++){
			double strain_gp[3];
			memset(strain_gp, 0, 3*sizeof(double));
			for(uint32_t k=0; k < elemtype->N_nodes; k++){
				strain_gp[0] += elemtype->Ni[k](elemtype->psi[j], elemtype->eta[j]) * 
					strain_on_nodes[connectivity_mtx[i*3+k] * 3];
				strain_gp[1] += elemtype->Ni[k](elemtype->psi[j], elemtype->eta[j]) * 
					strain_on_nodes[connectivity_mtx[i*3+k]*3+1];
				strain_gp[2] += elemtype->Ni[k](elemtype->psi[j], elemtype->eta[j]) * 
					strain_on_nodes[connectivity_mtx[i*3+k]*3+2];
			}
			uint32_t idx = i * N_gp + j;
			error[idx] += 
				sqrt(vcn_math_pow2(strain[idx * 3] - strain_gp[0]) +
				     vcn_math_pow2(strain[idx*3+1] - strain_gp[1]) +
				     vcn_math_pow2(strain[idx*3+2] - strain_gp[2]));
		}
	}

	/* Free memory */
	free(strain_on_nodes);
}

int vcn_fem_compute_2D_Solid_Mechanics
(const vcn_msh3trg_t *const mesh,
 const vcn_fem_elem_t *const elemtype,
 const vcn_fem_material_t *const material,
 const vcn_bcond_t *const bmeshcond,
 char enable_self_weight,
 double gravity[2],
 int (*solver)(const vcn_sparse_t *const A,
	       const double *const b,
	       double* x, uint32_t omp_threads),
 bool enable_plane_stress,
 double thickness,
 uint32_t omp_parallel_threads,
 bool* elements_enabled, /* NULL to enable all */
 double * displacement, /* Output */
 double * strain,       /* Output */
 const char* logfile    /* NULL: Does not create an output logifle */)
/* The output vectors must be allocated before start:
 *     > displacement:  2 * N_vertices (size of double)
 *     >       strain:  3 * N_vertices (size of double)
 */
{
	if(logfile != NULL){
		FILE *log = fopen(logfile, "a");
		fprintf(log, "Bidimensional elastic Solid Mechanics Solver using FEM\n");
		fclose(log);
	}

	/****************************************************************************/
	/********************** 1) Assemble system **********************************/
	/****************************************************************************/
	/* Allocate global Stiffness Matrix and Force Vector */
	vcn_graph_t *graph = vcn_msh3trg_create_vtx_graph(mesh);
	vcn_sparse_t *K = vcn_sparse_create(graph, NULL, 2);
	vcn_graph_destroy(graph);

	double* F = (double*)calloc(2 * mesh->N_vertices, sizeof(double));
	/* Allocate elemental Stiffness Matrix and Force Vector */
	pipeline_assemble_system(K, NULL, F,
				 mesh,
				 elemtype,
				 material,
				 enable_self_weight,
				 gravity,
				 enable_plane_stress,
				 thickness,
				 false, NULL,
				 elements_enabled);
	/*****************************************************************************/
	/********************** 2) Set boundary conditions ***************************/
	/*****************************************************************************/
	pipeline_set_boundary_conditions(K, F, bmeshcond, thickness, 1.0);

	/*****************************************************************************/
	/**************** 3) Solve system (to compute displacements) *****************/
	/*****************************************************************************/
  
	char solver_status = solver(K, F, displacement, omp_parallel_threads);
    
	/* Display failure info in logfile */
	if(solver_status != 0){
		if(logfile != NULL){
			FILE* log = fopen(logfile, "a");
			fprintf(log, "Solver fails (Code: %i).\n", solver_status);
			fclose(log);
		}
		vcn_sparse_destroy(K);
		free(F);
		return 1;
	}

	vcn_sparse_destroy(K);
	free(F);
	/*****************************************************************************/
	/********************* 4) Compute Strain            **************************/
	/*********************       >  S = B u             **************************/
	/*****************************************************************************/
	pipeline_compute_strain(strain, mesh, displacement, elemtype, false,
				enable_plane_stress, NULL, NULL, NULL, NULL);
	/* Successful exit */
	return 0;
}


void vcn_fem_compute_2D_Non_Linear_Solid_Mechanics
(const vcn_msh3trg_t *const mesh,
 const vcn_fem_elem_t *const elemtype,
 const vcn_fem_material_t *const material,
 const vcn_bcond_t *const bmeshcond,
 bool enable_self_weight,
 double gravity[2],
 bool enable_Cholesky_solver,
 bool enable_plane_stress,
 double thickness,
 vcn_fem_implicit_t* params,
 bool restore_computation,
 vcn_fem_output_t* output,
 const char* logfile)
/* Quasistatic formulation */
{
	uint32_t N_vertices = mesh->N_vertices;
	uint32_t N_elements = mesh->N_triangles;
	double* vertices = mesh->vertices;
	uint32_t* connectivity_mtx = mesh->vertices_forming_triangles;

	uint32_t omp_parallel_threads = 1;

	FILE *log = fopen(logfile, "a");
	fprintf(log, "FEM: Damage Model\n");
	fclose(log);
  
	uint32_t N_system_size = N_vertices * 2;

	uint32_t N_gp = elemtype->N_Gauss_points;

	/****************************************************************************/
	/*********************** > Allocate output **********************************/
	/****************************************************************************/
	double* displacement = (double*)calloc(2 * N_vertices, sizeof(double));
	double* strain = (double*)calloc(3 * N_elements * N_gp, sizeof(double));
	double* damage = (double*)calloc(N_elements * N_gp, sizeof(double));
	double* r_dmg = (double*)malloc(N_gp * N_elements * sizeof(double));

	/* Initialize r parameter used for damage calculation */
	for (uint32_t i = 0; i < N_gp * N_elements; i++)
		r_dmg[i] = material->r0_damage(material);

	/****************************************************************************/
	/*********************** > Allocate system **********************************/
	/****************************************************************************/
	/* Allocate global Stiffness Matrices */
	vcn_graph_t *graph = vcn_msh3trg_create_vtx_graph(mesh);
	vcn_sparse_t* K = vcn_sparse_create(graph, NULL, 2);
	vcn_sparse_t *L = NULL;
	vcn_sparse_t *U = NULL;
	/* Allocate the triangular matrices L and U using symbolic Cholesky */
	vcn_sparse_alloc_LU(K, &L, &U);

  
	/* Allocate force vectors and displacement increment */
	double* F = (double*)calloc(N_system_size, sizeof(double));
	double* P = (double*)calloc(N_system_size, sizeof(double));
	double* residual = (double*)calloc(N_system_size, sizeof(double));
	double* du = (double*)calloc(N_system_size, sizeof(double));

	/* Allocate damage parameter 'r' */
	double* r_dmg_prev = (double*)malloc(N_gp * N_elements * sizeof(double));
  
	/************************************************************************/
	/************************ > Start simulation of N steps *****************/
	/************************************************************************/
	for (uint32_t n = 0; n < vcn_fem_implicit_get_N_steps(params); n++) {
		log = fopen(logfile, "a");
		fprintf(log, "  [ Load step %i]\n", n + 1);
		fclose(log);
		memcpy(r_dmg_prev, r_dmg, N_gp * N_elements * sizeof(double));

		/***********************************************************************/
		/******************** > Implicit integration ***************************/
		/***********************************************************************/
		/* Implicit integration scheme */
		if(!enable_Cholesky_solver)
			memset(du, 0, N_system_size * sizeof(double));

		double residual_norm = 1;
		uint32_t residual_iter = 0;
		uint32_t residual_iter_without_enhance = 0;
		double residual_best;
		while(1) {
			/**********************************************************************/
			/******************** > Assemble system *******************************/
			/**********************************************************************/
			pipeline_assemble_system(K, NULL, F,
						 mesh,
						 elemtype,
						 material,
						 enable_self_weight,
						 gravity,
						 enable_plane_stress,
						 thickness,
						 true, /* Enable computing damage */
						 damage,
						 NULL);

			/**********************************************************************/
			/******************** > Set boundary conditions ***********************/
			/**********************************************************************/
			double condition_factor =
				(n + 1.0)/(double) vcn_fem_implicit_get_N_steps(params);

			/* Set Boundary Conditions */
			pipeline_set_boundary_conditions(K, F, bmeshcond, thickness, condition_factor);

			/**********************************************************************/
			/******************** > Verify residual *******************************/
			/**********************************************************************/
			/* Compute P increment */
			vcn_sparse_multiply_vector(K, displacement, P, omp_parallel_threads);

			/* Compute residual norm */
			residual_norm = 0;
			for(uint32_t i=0; i < N_system_size; i++){
				residual[i] = F[i] - P[i];
				residual_norm += vcn_math_pow2(residual[i]);
			}
      
			residual_norm = sqrt(residual_norm);

			/* Check if residual is minimized */
			if(residual_iter == 0){
				residual_best = residual_norm;
			}else{
				if(residual_norm < residual_best){
					residual_best = residual_norm;
					residual_iter_without_enhance = 0;
				}else
					residual_iter_without_enhance += 1;
			}

			/* Check stop criteria */
			if((residual_norm < params->residual_tolerance
			    && residual_iter > 0) || residual_iter >= params->N_max_iter ||
			   residual_iter_without_enhance >= params->N_max_iter_without_enhance){
				/* Write logfile */
				log = fopen(logfile, "a");
				fprintf(log, "    [%i: Residual: %e ]\n",
					residual_iter, residual_norm);
				fclose(log);
				break;
			}

			/**********************************************************************/
			/************** > Solve system (to compute displacements) *************/
			/**********************************************************************/
			char solver_status = 0; /* Default initialization (Must be initialized) */
			if(enable_Cholesky_solver){
				/* Decompose matrix */
				solver_status = 
					vcn_sparse_decompose_Cholesky(K, L, U, omp_parallel_threads);

				if(solver_status != 0)
					vcn_sparse_decompose_LU(K, L, U, omp_parallel_threads);
  
				/* Solve system */
				vcn_sparse_solve_LU(L, U, residual, du);
			}
			if(!enable_Cholesky_solver || solver_status != 0){
				if(solver_status != 0){
					log = fopen(logfile, "a");
					fprintf(log, "    [Cholesky Fails (The stiffness matrix is not positive definite).]\n");
					fprintf(log, "    [Solving with Conjugate Gradient (Tolerance: 1e-8)...]\n");
					fclose(log);
				}
				double error = 0.0; /* Default initialization (Must be initialized) */
				uint32_t iters = 0;     /* Default initialization (Must be initialized) */
				solver_status = 
					vcn_sparse_solve_CG_precond_Jacobi(K, residual, du,
									   vcn_sparse_get_size(K),
									   1e-8,
									   &iters, &error,
									   omp_parallel_threads);
				if(solver_status != 0){
					solver_status = vcn_sparse_solve_Gauss_Seidel(K, residual, du,
										      60 * vcn_sparse_get_size(K),
										      1e-8,
										      &iters, &error,
										      omp_parallel_threads);
					if(solver_status != 0){
						log = fopen(logfile, "a");
						fprintf(log, "    [Conjugate Gradient (p. Jacobi) and Gauss-Seidel fails.]\n");
						fprintf(log, "      [Gauss-Seidel iterations: %u]\n", iters);
						fprintf(log, "      [Gauss-Seidel error: %e ]\n", error);
						fclose(log);
					}
				}
			}
			/* Write logfile */
			log = fopen(logfile, "a");
			fprintf(log, "    [%i: Residual: %e]\n", residual_iter, residual_norm);
			fclose(log);

			/* Update displacements, velocities and accelerations */
			for(uint32_t i=0; i < N_system_size; i++)
				displacement[i] += du[i];
            

			/**********************************************************************/
			/**************** > Compute Strain and Damage         *****************/
			/**********************************************************************/
			pipeline_compute_strain(strain, mesh, displacement, elemtype, true,
						enable_plane_stress, material, damage,
						r_dmg_prev, r_dmg);

			/* Increase iterator */
			residual_iter ++;
		}
	}
	/*****************************************************************************/
	/***************************** > Free memory *********************************/
	/*****************************************************************************/
	free(displacement);
	free(strain);
	free(damage);

	vcn_sparse_destroy(K);
	if(enable_Cholesky_solver){
		vcn_sparse_destroy(L);
		vcn_sparse_destroy(U);
	}

	free(F);
	free(du);
	free(residual);
	free(P); 
  
	free(r_dmg_prev);
	free(r_dmg);

	if(restore_computation){
		free(vertices);
		free(connectivity_mtx);
	}
}

void vcn_fem_compute_stress_from_strain
(uint32_t N_elements,
 uint32_t* elements_connectivity_matrix, 
 const vcn_fem_elem_t *const elemtype,
 const vcn_fem_material_t *const material,
 bool enable_plane_stress,
 double* strain,
 uint32_t omp_parallel_threads,
 bool* elements_enabled /* NULL to enable all */,
 double* stress /* Output */){
	/* Compute stress from element strain */
#pragma omp parallel for num_threads(omp_parallel_threads) schedule(guided)
	for(uint32_t i = 0; i < N_elements; i++){
		/* Get material properties */
		double E = vcn_fem_material_get_elasticity_module(material);

		if(elements_enabled != NULL)
			if(!elements_enabled[i])
				E = 0;

		double v = vcn_fem_material_get_poisson_module(material);
    
		/* Get constitutive matrix */    
		double d11 = E/(1.0 - vcn_math_pow2(v));
		double d12 = v*d11;
    
		if(!enable_plane_stress){
			d11 = (E*(1.0-v))/((1.0 + v)*(1.0-2*v));
			d12 = (v*d11)/(1.0-v);
		}    
		double d22 = d11;
		double d33 = E/(2.0*(1.0+v));
		/* Calculate stress */
		stress[i * 3] = strain[i * 3] * d11 + strain[i*3+1] * d12;
		stress[i*3+1] = strain[i * 3] * d12 + strain[i*3+1] * d22;
		stress[i*3+2] = strain[i*3+2] * d33;
	}
}


void vcn_fem_interpolate_from_Gauss_points_to_nodes
(const vcn_msh3trg_t *const mesh, 
 const vcn_fem_elem_t *const elemtype,
 uint32_t N_components,
 double* values_on_GP_from_elements,
 double* values_interpolated_on_nodes /* Output */)
/* Let j be the elements adjacent to some vertex i.
 * Assuming dj is the distance from the vertex i to
 * the closest Gauss point of the element j.
 * pj is the value in such Gauss point and pi is the
 * interpolated value.
 *
 * Interpolation functions:
 *        ___
 *        \
 *  pi =  /__ wj pj  ,  wj = hj/bi, hj = 1/dj 
 *         j
 *
 *  where wj is the weight of the element j and
 *        ___
 *        \   1 /
 *  bi =  /__  / dk   <---- Normalization term
 *         k
 */
{
	uint32_t N_vertices = mesh->N_vertices;
	double *vertices = mesh->vertices;
	uint32_t N_elements = mesh->N_triangles;
	uint32_t* connectivity_mtx = mesh->vertices_forming_triangles;

	/* Define max number of elements for a single node */
	uint32_t N_room = 10;
	/* Allocate structures to store connectivities of nodes */
	uint32_t* N_elems_adjacents_to_node =
		(uint32_t*) calloc(N_vertices, sizeof(uint32_t));
	uint32_t* elems_adjacents_to_node =
		(uint32_t*) malloc(N_room * N_vertices * sizeof(uint32_t));
	/* Iterate over elements searching for connectivities */
	for(uint32_t i=0; i < N_elements; i++){
		/* Iterate over the nodes of the element */
		for(uint32_t j=0; j < elemtype->N_nodes; j++){
			uint32_t vj = connectivity_mtx[i*elemtype->N_nodes + j];
			elems_adjacents_to_node[vj*N_room + N_elems_adjacents_to_node[vj]] = i;
			if(N_elems_adjacents_to_node[vj] < N_room - 1)
				N_elems_adjacents_to_node[vj] += 1;
		}
	}
	/* Iterate over vertices to interpolate from the Elemental GP */
	memset(values_interpolated_on_nodes, 0,
	       N_components * N_vertices * sizeof(double));
	for(uint32_t i=0; i < N_vertices; i++){
		/* Iterate over the elements connected to the vertex */
		double sum_w = 0;
		for(uint32_t k=0; k < N_elems_adjacents_to_node[i]; k++){
			uint32_t elem_id = elems_adjacents_to_node[i*N_room + k];
			/* Get ID of the vertex relative to the element */
			uint32_t inside_idx = 0;
			while(connectivity_mtx[elem_id*elemtype->N_nodes + inside_idx])
				inside_idx++;    
			/* Get id of the closest GP */
			uint32_t id_gp =
				vcn_fem_elem_get_closest_Gauss_Point_to_the_ith_node(elemtype, inside_idx);

			/* Get coordinates to the GP */
			double gp[2] = {0,0};
			for(uint32_t j=0; j < elemtype->N_nodes; j++){
				uint32_t vj = connectivity_mtx[i*elemtype->N_nodes + j];
				gp[0] += vertices[vj * 2] *
					elemtype->Ni[j](elemtype->psi[id_gp], elemtype->eta[id_gp]);
				gp[1] += vertices[vj*2+1] *
					elemtype->Ni[j](elemtype->psi[id_gp], elemtype->eta[id_gp]);
			}
			/* Compute weight */
			double wk = 1.0 / vcn_utils2D_get_dist(gp, &(vertices[i*2]));
			sum_w += wk;
			/* Interpolate component by component */
			for(uint32_t c=0; c < N_components; c++){
				uint32_t cid = elem_id * elemtype->N_Gauss_points + id_gp;
				values_interpolated_on_nodes[i * N_components + c] += 
					wk * values_on_GP_from_elements[cid * N_components + c];
			}
      
		}
		/* Normalize component by component */
		for(uint32_t c=0; c < N_components; c++)
			values_interpolated_on_nodes[i * N_components + c] /= sum_w;
	}

	/* Free memory */
	free(N_elems_adjacents_to_node);
	free(elems_adjacents_to_node);
}
