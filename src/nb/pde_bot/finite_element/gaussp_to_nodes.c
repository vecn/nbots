#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <alloca.h>

#include "nb/memory_bot.h"
#include "nb/geometric_bot.h"
#include "nb/pde_bot/finite_element/element.h"
#include "nb/pde_bot/finite_element/gaussp_to_nodes.h"

#include "utils.h"

#define POW2(a) ((a)*(a))

static int assemble_system(const vcn_msh3trg_t *const mesh, 
			   const vcn_fem_elem_t *const elem,
			   uint32_t N_comp,
			   const double* gp_values,
			   double *M, double *b);

static int integrate_elemental_system
		       	(const vcn_fem_elem_t *elem,
			 uint32_t N_comp, const double *gp_values,
			 uint32_t id,
			 const vcn_msh3trg_t *mesh,
			 double *Me, double *be);

static void sum_gauss_point(const vcn_fem_elem_t *elem,
			    uint32_t N_comp, const double *gp_values,
			    uint32_t elem_id, int gp_id,
			    double detJ,
			    double *dNi_dx, double *dNi_dy,
			    double *Me, double *be);

static void add_to_global_system(const vcn_fem_elem_t *elem, 
				 uint32_t N_comp, uint32_t id,
				 const vcn_msh3trg_t *mesh,
				 double *Me, double *be,
				 double *M, double *b);

static void solve_system(const double *M, const double *b,
			 double *nodal_values,
			 uint32_t N_vtx, uint32_t N_comp);

static double get_gp_error(uint32_t id_elem, int id_gp,
			   uint32_t *conn_mtx,
			   const vcn_fem_elem_t *const elem,
			   uint32_t N_comp,
			   double* gp_values,
			   double* nodal_values);

int vcn_fem_interpolate_from_gpoints_to_nodes
		(const vcn_msh3trg_t *const mesh, 
		 const vcn_fem_elem_t *const elem,
		 uint32_t N_comp,
		 const double* gp_values,
		 double* nodal_values /* Output */)
{
	int status = 1;
	double* M = calloc(mesh->N_vertices, sizeof(*M));
	double* b = calloc(N_comp * mesh->N_vertices, sizeof(*b));
	status = assemble_system(mesh, elem, N_comp,
				 gp_values, M, b);

	if (0 != status)
		goto CLEANUP;

	solve_system(M, b, nodal_values, mesh->N_vertices, N_comp);

	status = 0;
CLEANUP:
	free(M);
	free(b);
	return status;
}

static int assemble_system(const vcn_msh3trg_t *const mesh, 
			   const vcn_fem_elem_t *const elem,
			   uint32_t N_comp,
			   const double* gp_values,
			   double *M, double *b)
{
	int status = 1;

	uint32_t N_elem = mesh->N_triangles;

	uint8_t N_nodes_x_elem = vcn_fem_elem_get_N_nodes(elem);
	
	uint16_t nodes_size = N_nodes_x_elem * sizeof(double);
	uint32_t memsize = (1 + N_comp) * nodes_size;
	char *memblock = NB_SOFT_MALLOC(memsize);
	double* Me = (void*) memblock;
	double* be = (void*) (memblock + nodes_size);

	for (uint32_t i = 0; i < N_elem; i++) {
		status = integrate_elemental_system(elem, N_comp,
						    gp_values,
						    i, mesh,
						    Me, be);

		if (0 != status)
			goto EXIT;

		add_to_global_system(elem, N_comp, i, mesh,
				     Me, be, M, b);
	}
EXIT:
	NB_SOFT_FREE(memsize, memblock);
	return status;
}

static int integrate_elemental_system
		       	(const vcn_fem_elem_t *elem,
			 uint32_t N_comp, const double *gp_values,
			 uint32_t id,
			 const vcn_msh3trg_t *mesh,
			 double *Me, double *be)
{
	int status = 1;

	uint8_t N_nodes = vcn_fem_elem_get_N_nodes(elem);
	memset(Me, 0, N_nodes * sizeof(*Me));
	memset(be, 0, N_comp * N_nodes * sizeof(*be));

	/* Allocate Cartesian derivatives for each Gauss Point */
	uint32_t deriv_memsize = N_nodes * sizeof(double);
	char *deriv_memblock = NB_SOFT_MALLOC(2 * deriv_memsize);
	double *dNi_dx = (void*) (deriv_memblock);
	double *dNi_dy = (void*) (deriv_memblock + deriv_memsize);

	uint8_t N_gp = vcn_fem_elem_get_N_gpoints(elem);
	for (uint32_t j = 0; j < N_gp; j++) {
		double Jinv[4];
		double detJ = nb_fem_get_jacobian(elem, id,
						  mesh, j, Jinv);

		if (nb_fem_elem_is_distorted(detJ))
			goto EXIT;

		nb_fem_get_derivatives(elem, j, Jinv, dNi_dx, dNi_dy);

		sum_gauss_point(elem, N_comp, gp_values, id, j, detJ,
				dNi_dx, dNi_dy, Me, be);
	}
	status = 0;
EXIT:
	NB_SOFT_FREE(2 * deriv_memsize, deriv_memblock);
	return status;
}

static void sum_gauss_point(const vcn_fem_elem_t *elem,
			    uint32_t N_comp, const double *gp_values,
			    uint32_t elem_id, int gp_id,
			    double detJ,
			    double *dNi_dx, double *dNi_dy,
			    double *Me, double *be)
{
	double wp = vcn_fem_elem_weight_gp(elem, gp_id);
	uint8_t N_nodes = vcn_fem_elem_get_N_nodes(elem);
	for (uint8_t i = 0; i < N_nodes; i++) {
		double Ni = vcn_fem_elem_Ni(elem, i, gp_id);
		for (uint32_t j = 0; j < N_nodes; j++) {
			double Nj = vcn_fem_elem_Ni(elem, j, gp_id);
			double integral = Ni * Nj * detJ * wp;
			Me[i] += integral;
		}

		uint8_t N_gp = vcn_fem_elem_get_N_gpoints(elem);
		uint32_t global_gp = elem_id * N_gp + gp_id;
		double integral = Ni * detJ * wp;
		for (int c = 0; c < N_comp; c++) {
			double val = gp_values[global_gp * N_comp + c];
			be[i * N_comp + c] += val * integral;
		}
	}
}

static void add_to_global_system(const vcn_fem_elem_t *elem, 
				 uint32_t N_comp, uint32_t id,
				 const vcn_msh3trg_t *mesh,
				 double *Me, double *be,
				 double *M, double *b)
{
	uint8_t N_nodes = vcn_fem_elem_get_N_nodes(elem);
	for (uint32_t i = 0; i < N_nodes; i++) {
		uint32_t v1 = mesh->vertices_forming_triangles[id*3+i];
		M[v1] += Me[i];

		for (int c = 0; c < N_comp; c++)
			b[v1 * N_comp + c] += be[i * N_comp + c];
	}
}

static void solve_system(const double *M, const double *b,
			 double *nodal_values,
			 uint32_t N_vtx, uint32_t N_comp)
{
	for (uint32_t i = 0; i < N_vtx; i++) {
		for (uint32_t c = 0; c < N_comp; c++) {
			uint32_t id = i * N_comp + c;
			nodal_values[id] = b[id] / M[i];
		}
	}
}

void vcn_fem_get_error_on_gpoints(const vcn_msh3trg_t *const mesh,
				  const vcn_fem_elem_t *const elem,
				  uint32_t N_comp,
				  double* gp_values,
				  double* gp_error)
{
	uint32_t N_vtx = mesh->N_vertices;
	uint32_t N_elem = mesh->N_triangles;

	uint32_t N_gp = vcn_fem_elem_get_N_gpoints(elem);
	memset(gp_error, 0, N_elem * N_gp * sizeof(double));

	/* Interpolate strain to nodes */
	double* nodal_values = calloc(N_comp * N_vtx, sizeof(double));
	vcn_fem_interpolate_from_gpoints_to_nodes(mesh, elem,
						  N_comp, gp_values,
						  nodal_values);

	uint32_t *conn_mtx = mesh->vertices_forming_triangles;

	for (uint32_t i = 0; i < N_elem; i++) {
		for (uint32_t j = 0; j < N_gp; j++) {
			gp_error[i * N_gp + j] = get_gp_error(i, j, conn_mtx,
							      elem, N_comp,
							      gp_values,
							      nodal_values);
		}
	}

	/* Free memory */
	free(nodal_values);
}

static double get_gp_error(uint32_t id_elem, int id_gp,
			   uint32_t *conn_mtx,
			   const vcn_fem_elem_t *const elem,
			   uint32_t N_comp,
			   double* gp_values,
			   double* nodal_values)
{
	double *value_gp = alloca(N_comp * sizeof(*value_gp));
	memset(value_gp, 0, N_comp * sizeof(*value_gp));

	uint8_t N_nodes = vcn_fem_elem_get_N_nodes(elem);
	for (uint32_t i = 0; i < N_nodes; i++) {
		uint32_t vid = conn_mtx[id_elem * N_nodes + i];
		double Ni = vcn_fem_elem_Ni(elem, i, id_gp);
		for (int c = 0; c < N_comp; c++)
			value_gp[c] += Ni * nodal_values[vid * N_comp + c];
	}
	uint8_t N_gp = vcn_fem_elem_get_N_gpoints(elem);
	uint32_t idx = id_elem * N_gp + id_gp;
	double sum = 0.0;
	for (int c = 0; c < N_comp; c++) {
		double real_value = gp_values[idx * N_comp + c];
		double error;
		if (fabs(real_value) > 1e-10) 
			error = fabs(1.0 - value_gp[c]/real_value);
		else
			error = fabs(value_gp[c]);
		sum += POW2(error);
	}
	return sqrt(sum);
}
