#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <alloca.h>

#include "nb/math_bot.h"
#include "nb/cfreader_cat.h"
#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/eigen_bot.h"
#include "nb/geometric_bot.h"
#include "nb/graph_bot.h"
#include "nb/pde_bot/material.h"
#include "nb/pde_bot/boundary_conditions/bcond.h"
#include "nb/pde_bot/boundary_conditions/bcond_iter.h"
#include "nb/pde_bot/finite_element/element.h"
#include "nb/pde_bot/finite_element/gaussp_to_nodes.h"
#include "nb/pde_bot/finite_element/solid_mechanics/static_elasticity2D.h"

#include "../utils.h"
#include "pipeline.h"

#define POW2(a) ((a)*(a))

static int assemble_element(const vcn_fem_elem_t *elem, uint32_t id,
			    const vcn_msh3trg_t *mesh,
			    const vcn_fem_material_t *material,
			    bool is_enabled,
			    nb_analysis2D_t analysis2D,
			    nb_analysis2D_params *params2D,
			    bool enable_self_weight,
			    double gravity[2],
			    vcn_sparse_t *K, double *M, double *F);
static void set_plane_stress(double D[4], double E, double v);
static void set_plane_strain(double D[4], double E, double v);
static int integrate_elemental_system
		       	(const vcn_fem_elem_t *elem, uint32_t id,
			 double D[4], double density, double gravity[2],
			 const vcn_msh3trg_t *mesh,
			 nb_analysis2D_params *params2D,
			 bool enable_self_weight,
			 double *Ke, double *Me, double *Fe);
static void set_bcond_neumann_sgm(const vcn_msh3trg_t *msh3trg,
				  vcn_sparse_t* K, double* F, 
				  const nb_bcond_t *const bcond, 
				  double factor);

static double get_input_sgm_length(const vcn_msh3trg_t *msh3trg,
				   uint32_t sgm_id);
static double get_input_subsgm_length(const vcn_msh3trg_t *msh3trg,
				      uint32_t sgm_id, uint32_t subsgm_id);
static void set_bcond_neumann(const double *vertices,
			      uint8_t N_dof,
			      double* F, double factor,
			      nb_bcond_iter_t *iter,
			      uint32_t vtx_id);
static void set_bcond_neumann_vtx(const vcn_msh3trg_t *msh3trg,
				  vcn_sparse_t* K, double* F, 
				  const nb_bcond_t *const bcond, 
				  double factor);
static void set_bcond_dirichlet_sgm(const vcn_msh3trg_t *msh3trg,
				    vcn_sparse_t* K, double* F, 
				    const nb_bcond_t *const bcond, 
				    double factor);
static void set_bcond_dirichlet(const double *vertices,
				vcn_sparse_t* K, uint8_t N_dof,
				double* F, double factor,
				nb_bcond_iter_t *iter,
				uint32_t vtx_id);
static void set_bcond_dirichlet_vtx(const vcn_msh3trg_t *msh3trg,
				    vcn_sparse_t* K, double* F, 
				    const nb_bcond_t *const bcond, 
				    double factor);
static int get_element_strain(uint32_t id, double *strain,
			      const vcn_msh3trg_t *const mesh,
			      double *displacement,
			      const vcn_fem_elem_t *const elem,
			      nb_analysis2D_t analysis2D,
			      const vcn_fem_material_t *const material);

int pipeline_assemble_system
		(vcn_sparse_t* K, double* M, double *F,
		 const vcn_msh3trg_t *const mesh,
		 const vcn_fem_elem_t *const elem,
		 const vcn_fem_material_t *const material,
		 bool enable_self_weight,
		 double gravity[2],
		 nb_analysis2D_t analysis2D,
		 nb_analysis2D_params *params2D,
		 const bool* elements_enabled /* NULL to enable all */)
{
	int status = 1;
	uint32_t N_elem = mesh->N_triangles;
	vcn_sparse_reset(K);
	if (NULL != M)
		memset(M, 0, vcn_sparse_get_size(K) * sizeof(*M));
	memset(F, 0, vcn_sparse_get_size(K) * sizeof(*F));

	for (uint32_t i = 0; i < N_elem; i++) {
		bool is_enabled = pipeline_elem_is_enabled(elements_enabled, i);
		int status_element =
			assemble_element(elem, i, mesh, material, is_enabled,
					 analysis2D, params2D,
					 enable_self_weight, gravity,
					 K, M, F);
		if (0 != status_element)
			goto EXIT;
	}
	status = 0;
EXIT:
	return status;
}

bool pipeline_elem_is_enabled(const bool *elements_enabled, uint32_t id)
{
	bool is_enabled = true;
	if (NULL != elements_enabled)
		is_enabled = elements_enabled[id];
	return is_enabled;
}

static int assemble_element(const vcn_fem_elem_t *elem, uint32_t id,
			    const vcn_msh3trg_t *mesh,
			    const vcn_fem_material_t *material,
			    bool is_enabled,
			    nb_analysis2D_t analysis2D,
			    nb_analysis2D_params *params2D,
			    bool enable_self_weight,
			    double gravity[2],
			    vcn_sparse_t *K, double *M, double *F)
{
	double D[4] = {1e-6, 1e-6, 1e-6, 1e-6};
	double density = 1e-6;
	if (is_enabled) {
		pipeline_get_constitutive_matrix(D, material, analysis2D);
		density = vcn_fem_material_get_density(material);
	}

	uint8_t N_nodes = vcn_fem_elem_get_N_nodes(elem);
	double* Ke = malloc(4 * POW2(N_nodes) * sizeof(*Ke));
	double* Me = NULL;
	if(M != NULL)
		Me = malloc(2 * N_nodes * sizeof(*Me));
	double* Fe = malloc(2 * N_nodes * sizeof(*Fe));
	
	int status = integrate_elemental_system(elem, id, D, density,
						gravity, mesh, params2D,
						enable_self_weight,
						Ke, Me, Fe);
	if (0 != status)
		goto CLEANUP;
	
	pipeline_add_to_global_system(elem, id, mesh, Ke, Me, Fe, K, M, F);

CLEANUP:
	free(Ke);
	if (NULL != M) 
		free(Me);
	free(Fe);
	return status;
}

void pipeline_get_constitutive_matrix(double D[4], 
				      const vcn_fem_material_t *material,
				      nb_analysis2D_t analysis2D)
{
	double E = vcn_fem_material_get_elasticity_module(material);
	double v = vcn_fem_material_get_poisson_module(material);
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

static int integrate_elemental_system
		       	(const vcn_fem_elem_t *elem, uint32_t id,
			 double D[4], double density, double gravity[2],
			 const vcn_msh3trg_t *mesh,
			 nb_analysis2D_params *params2D,
			 bool enable_self_weight,
			 double *Ke, double *Me, double *Fe)
{
	int status = 1;

	double fx = 0.0;
	double fy = 0.0;
	if (enable_self_weight) {
		fx = gravity[0] * density;
		fy = gravity[1] * density;
	}
    
	uint8_t N_nodes = vcn_fem_elem_get_N_nodes(elem);
	memset(Ke, 0, 4 * POW2(N_nodes) * sizeof(*Ke));
	if(Me != NULL)
		memset(Me, 0, 2 * N_nodes * sizeof(*Me));
	memset(Fe, 0, 2 * N_nodes * sizeof(*Fe));

	/* Allocate Cartesian derivatives for each Gauss Point */
	uint32_t deriv_memsize = N_nodes * sizeof(double);
	char *deriv_memblock = NB_SOFT_MALLOC(2 * deriv_memsize);
	double *dNi_dx = (void*) (deriv_memblock);
	double *dNi_dy = (void*) (deriv_memblock + deriv_memsize);

	uint8_t N_gp = vcn_fem_elem_get_N_gpoints(elem);
	for (uint32_t j = 0; j < N_gp; j++) {
		double Jinv[4];
		double detJ = nb_fem_get_jacobian(elem, id, mesh, j, Jinv);

		if (nb_fem_elem_is_distorted(detJ))
			goto EXIT;

		nb_fem_get_derivatives(elem, j, Jinv, dNi_dx, dNi_dy);

		double thickness = params2D->thickness;
		pipeline_sum_gauss_point(elem, j, D, density, thickness,
					 detJ, dNi_dx, dNi_dy, fx, fy,
					 Ke, Me, Fe);
	}
	status = 0;
EXIT:
	NB_SOFT_FREE(2 * deriv_memsize, deriv_memblock);
	return status;
}

void pipeline_sum_gauss_point(const vcn_fem_elem_t *elem, int gp_id,
			      double D[4], double density,
			      double thickness, double detJ,
			      double *dNi_dx, double *dNi_dy,
			      double fx, double fy,
			      double *Ke, double *Me, double *Fe)
{
	/* Compute elemental stiffness matrix
	 *       _            _         _        _        _  _
	 *      |dNi/dx        |       | D0 D1    |      |  |
	 * Bi = |       dNi/dy |   D = | D1 D2    |  K = |  | B'DB dx dy
	 *      |dNi/dy dNi/dx_|       |_      D3_|     _| _|

	 * B  = [B1 B2...Bi... Bn]
	 *
	 */
	double wp = vcn_fem_elem_weight_gp(elem, gp_id);
      
	uint8_t N_nodes = vcn_fem_elem_get_N_nodes(elem);
	for (uint32_t i = 0; i < N_nodes; i++) {
		double Ni = vcn_fem_elem_Ni(elem, i, gp_id);
		for (uint32_t j = 0; j < N_nodes; j++) {
			/*  Integrating elemental siffness matrix */
			Ke[(i * 2)*(2 * N_nodes) + (j * 2)] += 
				(dNi_dx[i]*dNi_dx[j]*D[0] +
				 dNi_dy[i]*dNi_dy[j]*D[3]) *
				detJ * thickness * wp;
			Ke[(i * 2)*(2 * N_nodes) + (j*2+1)] +=
				(dNi_dx[i]*dNi_dy[j]*D[1] +
				 dNi_dy[i]*dNi_dx[j]*D[3]) *
				detJ * thickness * wp;
			Ke[(i*2+1)*(2 * N_nodes) + (j * 2)] +=
				(dNi_dy[i]*dNi_dx[j]*D[1] +
				 dNi_dx[i]*dNi_dy[j]*D[3]) *
				detJ * thickness * wp;
			Ke[(i*2+1)*(2 * N_nodes) + (j*2+1)] +=
				(dNi_dy[i]*dNi_dy[j]*D[2] +
				 dNi_dx[i]*dNi_dx[j]*D[3]) *
				detJ * thickness * wp;

			/*  Integrating elemental mass matrix */
			double Nj = vcn_fem_elem_Ni(elem, j, gp_id);
			if (NULL != Me) {
				/* OPPORTUNITY: Allocate one for each comp */
				double integral = Ni * Nj * density * detJ *
					thickness * wp;
				Me[i * 2] += integral;
				Me[i*2+1] += integral;
			}
		}
		/* Compute elemental forces vector */
		/* OPPORTUNITY: Allocate just one for each component */
		double integral = Ni * detJ * thickness * wp;
		Fe[i * 2] += integral * fx;
		Fe[i*2+1] += integral * fy;
	}
}

void pipeline_add_to_global_system(const vcn_fem_elem_t *elem, uint32_t id,
				   const vcn_msh3trg_t *mesh,
				   double *Ke, double *Me, double *Fe,
				   vcn_sparse_t *K, double *M, double *F)
{
	uint8_t N_nodes = vcn_fem_elem_get_N_nodes(elem);
	for (uint32_t i = 0; i < N_nodes; i++) {
		uint32_t v1 = mesh->vertices_forming_triangles[id*3+i];
		for (uint32_t j = 0; j < N_nodes; j++) {
			uint32_t v2 = mesh->vertices_forming_triangles[id*3+j];
			vcn_sparse_add(K, v1 * 2, v2 * 2,
				       Ke[(i * 2)*(2 * N_nodes) +
					  (j * 2)]);
			vcn_sparse_add(K, v1*2 + 1, v2 * 2,
				       Ke[(i*2 + 1)*(2 * N_nodes) +
					  (j * 2)]);
			vcn_sparse_add(K, v1 * 2, v2*2 + 1,
				       Ke[(i * 2)*(2 * N_nodes) +
					  (j*2+1)]);
			vcn_sparse_add(K, v1*2 + 1, v2*2 + 1,
				       Ke[(i*2 + 1)*(2 * N_nodes) +
					  (j*2+1)]);
		}
		/* Add to global mass matrix */
		if (NULL != M) {
			M[v1 * 2] += Me[i * 2];
			M[v1*2+1] += Me[i*2+1];
		}
		/* Add to global forces vector */
		F[v1 * 2] += Fe[i * 2];
		F[v1*2+1] += Fe[i*2+1];
	}
}

void pipeline_set_boundary_conditions(const vcn_msh3trg_t *msh3trg,
				      vcn_sparse_t* K, double* F, 
				      const nb_bcond_t *const bcond,
				      double factor)
{
	set_bcond_neumann_sgm(msh3trg, K, F, bcond, factor);
	set_bcond_neumann_vtx(msh3trg, K, F, bcond, factor);
	set_bcond_dirichlet_sgm(msh3trg, K, F, bcond, factor);
	set_bcond_dirichlet_vtx(msh3trg, K, F, bcond, factor);
}

static void set_bcond_neumann_sgm(const vcn_msh3trg_t *msh3trg,
				  vcn_sparse_t* K,
				  double* F, 
				  const nb_bcond_t *const bcond, 
				  double factor)
{
	uint8_t N_dof = nb_bcond_get_N_dof(bcond);
	uint16_t size = nb_bcond_iter_get_memsize();
	nb_bcond_iter_t *iter = alloca(size);
	nb_bcond_iter_init(iter);
	nb_bcond_iter_set_conditions(iter, bcond, NB_NEUMANN,
				     NB_BC_ON_SEGMENT);
	while (nb_bcond_iter_has_more(iter)) {
		nb_bcond_iter_go_next(iter);
		uint32_t model_id = nb_bcond_iter_get_id(iter);
		double sgm_length = get_input_sgm_length(msh3trg, model_id);

		uint32_t N = msh3trg->N_vtx_x_inputsgm[model_id];
		for (uint32_t i = 0; i < N - 1; i++) {
			double subsgm_length =
				get_input_subsgm_length(msh3trg, model_id, i);
			double w = subsgm_length / sgm_length;

			uint32_t v1_id = 
				msh3trg->meshvtx_x_inputsgm[model_id][i];
			set_bcond_neumann(msh3trg->vertices, N_dof, F,
					  factor * w * 0.5, iter, v1_id);/* AQUI VOY */

			uint32_t v2_id = 
				msh3trg->meshvtx_x_inputsgm[model_id][i + 1];
			set_bcond_neumann(msh3trg->vertices, N_dof, F,
					  factor * w * 0.5, iter, v2_id);
		}
	}
	nb_bcond_iter_finish(iter);
}

static double get_input_sgm_length(const vcn_msh3trg_t *msh3trg,
				   uint32_t sgm_id)
{
	uint32_t last_vtx = msh3trg->N_vtx_x_inputsgm[sgm_id] - 1;
	uint32_t v1 = msh3trg->meshvtx_x_inputsgm[sgm_id][0];
	uint32_t v2 = msh3trg->meshvtx_x_inputsgm[sgm_id][last_vtx];
	double x1 = msh3trg->vertices[v1 * 2];
	double y1 = msh3trg->vertices[v1*2+1];
	double x2 = msh3trg->vertices[v2 * 2];
	double y2 = msh3trg->vertices[v2*2+1];
	return sqrt(POW2(x1 - x2) + POW2(y1 - y2));
}

static double get_input_subsgm_length(const vcn_msh3trg_t *msh3trg,
				      uint32_t sgm_id, uint32_t subsgm_id)
{
	uint32_t v1 = msh3trg->meshvtx_x_inputsgm[sgm_id][subsgm_id];
	uint32_t v2 = msh3trg->meshvtx_x_inputsgm[sgm_id][subsgm_id + 1];
	double x1 = msh3trg->vertices[v1 * 2];
	double y1 = msh3trg->vertices[v1*2+1];
	double x2 = msh3trg->vertices[v2 * 2];
	double y2 = msh3trg->vertices[v2*2+1];
	return sqrt(POW2(x1 - x2) + POW2(y1 - y2));
}

static void set_bcond_neumann(const double *vertices,
			      uint8_t N_dof,
			      double* F, double factor,
			      nb_bcond_iter_t *iter,
			      uint32_t vtx_id)
{
	for (uint8_t j = 0; j < N_dof; j++) {
		bool mask = nb_bcond_iter_get_mask(iter, j);
		double *x = &(vertices[vtx_id * 2]);
		double val = nb_bcond_iter_get_val(iter, j, x, 0);
		if (mask) {
			uint32_t mtx_id = vtx_id * N_dof + j;
			F[mtx_id] += factor * val;
		}
	}

}

static void set_bcond_neumann_vtx(const vcn_msh3trg_t *msh3trg,
				  vcn_sparse_t* K,
				  double* F, 
				  const nb_bcond_t *const bcond,
				  double factor)
{
	uint8_t N_dof = nb_bcond_get_N_dof(bcond);
	uint16_t size = nb_bcond_iter_get_memsize();
	nb_bcond_iter_t *iter = alloca(size);
	nb_bcond_iter_init(iter);
	nb_bcond_iter_set_conditions(iter, bcond, NB_NEUMANN,
				     NB_BC_ON_POINT);
	while (nb_bcond_iter_has_more(iter)) {
		nb_bcond_iter_go_next(iter);
		uint32_t model_id = nb_bcond_iter_get_id(iter);
		uint32_t mesh_id = msh3trg->input_vertices[model_id];
		set_bcond_neumann(msh3trg->vertices, N_dof, F,
				  factor, iter, mesh_id);
	}
	nb_bcond_iter_finish(iter);
}

static void set_bcond_dirichlet_sgm(const vcn_msh3trg_t *msh3trg,
				    vcn_sparse_t* K,
				    double* F, 
				    const nb_bcond_t *const bcond,
				    double factor)
{
	uint8_t N_dof = nb_bcond_get_N_dof(bcond);
	uint16_t size = nb_bcond_iter_get_memsize();
	nb_bcond_iter_t *iter = alloca(size);
	nb_bcond_iter_init(iter);
	nb_bcond_iter_set_conditions(iter, bcond, NB_DIRICHLET,
				     NB_BC_ON_SEGMENT);
	while (nb_bcond_iter_has_more(iter)) {
		nb_bcond_iter_go_next(iter);
		uint32_t model_id = nb_bcond_iter_get_id(iter);
		uint32_t N = msh3trg->N_vtx_x_inputsgm[model_id];
		for (uint32_t i = 0; i < N; i++) {
			uint32_t mesh_id = 
				msh3trg->meshvtx_x_inputsgm[model_id][i];
			set_bcond_dirichlet(msh3trg->vertices, 
					    K, N_dof, F, factor,
					    iter, mesh_id);
		}
	}
	nb_bcond_iter_finish(iter);
}

static void set_bcond_dirichlet(const double *vertices,
				vcn_sparse_t* K, uint8_t N_dof,
				double* F, double factor,
				nb_bcond_iter_t *iter,
				uint32_t vtx_id)
{
	for (uint8_t j = 0; j < N_dof; j++) {
		bool mask = nb_bcond_iter_get_mask(iter, j);
		double *x = &(vertices[vtx_id * 2]);
		double val = nb_bcond_iter_get_val(iter, j, x, 0);
		if (mask) {
			uint32_t mtx_id = vtx_id * N_dof + j;
			vcn_sparse_set_Dirichlet_condition(K, F, mtx_id,
							   factor * val);
		}
	}
}

static void set_bcond_dirichlet_vtx(const vcn_msh3trg_t *msh3trg,
				    vcn_sparse_t* K,
				    double* F, 
				    const nb_bcond_t *const bcond,
				    double factor)
{
	uint8_t N_dof = nb_bcond_get_N_dof(bcond);
	uint16_t size = nb_bcond_iter_get_memsize();
	nb_bcond_iter_t *iter = alloca(size);
	nb_bcond_iter_init(iter);
	nb_bcond_iter_set_conditions(iter, bcond, NB_DIRICHLET,
				     NB_BC_ON_POINT);
	while (nb_bcond_iter_has_more(iter)) {
		nb_bcond_iter_go_next(iter);
		uint32_t model_id = nb_bcond_iter_get_id(iter);
		uint32_t mesh_id = msh3trg->input_vertices[model_id];
		set_bcond_dirichlet(msh3trg->vertices, K, N_dof, F,
				    factor, iter, mesh_id);
	}
	nb_bcond_iter_finish(iter);
}

void pipeline_compute_strain(double *strain,
			     const vcn_msh3trg_t *const mesh,
			     double *displacement,
			     const vcn_fem_elem_t *const elem,
			     nb_analysis2D_t analysis2D,
			     const vcn_fem_material_t *const material)
{
	uint32_t N_elem = mesh->N_triangles;

	uint8_t N_gp = vcn_fem_elem_get_N_gpoints(elem);
	memset(strain, 0, 3 * N_gp * N_elem * sizeof(*strain));

	for (uint32_t i = 0; i < N_elem; i++)
		get_element_strain(i, strain, mesh, displacement,
				   elem, analysis2D, material);
}

static int get_element_strain(uint32_t id, double *strain,
			      const vcn_msh3trg_t *const mesh,
			      double *displacement,
			      const vcn_fem_elem_t *const elem,
			      nb_analysis2D_t analysis2D,
			      const vcn_fem_material_t *const material)
{
	int status = 1;
	uint32_t *conn_mtx = mesh->vertices_forming_triangles;

	/* Allocate Cartesian derivatives for each Gauss Point */
	uint8_t N_nodes = vcn_fem_elem_get_N_nodes(elem);
	uint32_t deriv_memsize = N_nodes * sizeof(double);
	char *deriv_memblock = NB_SOFT_MALLOC(2 * deriv_memsize);
	double *dNi_dx = (void*) (deriv_memblock);
	double *dNi_dy = (void*) (deriv_memblock + deriv_memsize);

	/* Integrate domain */
	uint8_t N_gp = vcn_fem_elem_get_N_gpoints(elem);
	for (uint32_t j = 0; j < N_gp; j++) {
		double Jinv[4];
		double detJ = nb_fem_get_jacobian(elem, id, mesh, j, Jinv);
      
		if (nb_fem_elem_is_distorted(detJ))
			goto EXIT;

		nb_fem_get_derivatives(elem, j, Jinv, dNi_dx, dNi_dy);

		uint32_t idx = id * N_gp + j;
		/* Compute Strain at Gauss Point */
		for (uint32_t i = 0; i < N_nodes; i++) {
			uint32_t inode = conn_mtx[id * N_nodes + i];
			strain[idx * 3] += dNi_dx[i] * displacement[inode * 2];
			strain[idx*3+1] += dNi_dy[i] * displacement[inode*2+1];
			strain[idx*3+2] += (dNi_dy[i] * displacement[inode * 2] +
					    dNi_dx[i] * displacement[inode*2+1]);
		}
	}
	status = 0;
EXIT:
	NB_SOFT_FREE(2 * deriv_memsize, deriv_memblock);
	return status;
}

void pipeline_compute_main_stress(double *stress, 
				  double *main_stress,
				  uint32_t N_elements,
				  const vcn_fem_elem_t *const elem)
{
	uint8_t N_gp = vcn_fem_elem_get_N_gpoints(elem);
	for(uint32_t i=0; i < N_elements; i++){
		for(uint32_t j=0; j < N_gp; j++){
			uint32_t idx = i * N_gp + j;
			double sigma_avg = 0.5 * (stress[idx * 3] +
						  stress[idx*3+1]);
			double R = sqrt(0.25 * 
					POW2(stress[idx * 3] -
					     stress[idx*3+1]) +
					POW2(stress[idx*3+2]));
			main_stress[idx * 2] = sigma_avg + R;
			main_stress[idx*2+1] = sigma_avg - R;
		}
	}
}
