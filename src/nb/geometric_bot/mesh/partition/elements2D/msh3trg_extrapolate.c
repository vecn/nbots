#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <alloca.h>

#include "nb/memory_bot.h"
#include "nb/geometric_bot/mesh/mshition/msh3trg.h"

#define POW2(a) ((a)*(a))

static void assemble_system(const void *msh, 
			    const double* elem_values,
			    uint32_t N_comp,
			    double *M, double *f);

static void assemble_elem(const void *msh,
			  const double *elem_values,
			  uint32_t elem_id,
			  uint8_t N_comp,
			  double *M, double *F);
static double get_jacobian(const void *msh, uint32_t id, double Jinv[4]);
static void get_deriv(double Jinv[4], double dNi_dx[3], double dNi_dy[3]);
static int integrate_elemental_system
		       	(const vcn_fem_elem_t *elem,
			 uint32_t N_comp, const double *gp_values,
			 uint32_t id,
			 const void *msh,
			 double *Me, double *be);

static void sum_gauss_point(const vcn_fem_elem_t *elem,
			    uint32_t N_comp, const double *gp_values,
			    uint32_t elem_id, int gp_id,
			    double detJ,
			    double *dNi_dx, double *dNi_dy,
			    double *Me, double *be);

static void solve_system(const double *M, const double *b,
			 double *nodal_values,
			 uint32_t N_vtx, uint32_t N_comp);

static double get_gp_error(uint32_t id_elem, int id_gp,
			   const void *msh,
			   const vcn_fem_elem_t *const elem,
			   uint32_t N_comp,
			   double* gp_values,
			   double* nodal_values);


void nb_msh3trg_extrapolate_elems_to_nodes(const void *msh, uint8_t N_comp,
					   const double *elem_values,
					   double *nodal_values)
{
	uint32_t N_nodes = nb_msh3trg_get_N_nodes(msh);

	char *memblock = malloc((1 + N_comp) * N_nodes * sizeof(double));

	double *M = (void*) memblock;
	double *F = (void*) (memblock + N_nodes * sizeof(double));
	memset(M, 0, N_nodes * sizeof(*M));
	memset(F, 0, N_nodes * N_comp * sizeof(*F));
	
	assemble_system(msh, elem_values, N_comp, M, F);

	solve_system(msh, N_comp, M, F, nodal_values);

	free(memblock);
}

static void assemble_system(const void *msh, 
			    const double* elem_values,
			    uint32_t N_comp,
			    double *M, double *f)
{
	uint32_t N_elems = nb_msh3trg_get_N_elems(msh);
	for (uint32_t i = 0; i < N_elems; i++)
		assemble_elem(msh, elem_values, i, N_comp, M, F);
}

static void assemble_elem(const void *msh,
			  const double *elem_values,
			  uint32_t elem_id,
			  uint8_t N_comp,
			  double *M, double *F)
{
	double Jinv[4];/* AQUI VOY, ya solo este y el de mshquad */
	double detJ = nb_fem_get_jacobian(elem, id, msh, 0, Jinv);

	double dNi_dx[3];
	double dNi_dy[3];
	nb_fem_get_derivatives(elem, j, Jinv, dNi_dx, dNi_dy);

	sum_gauss_point(elem, N_comp, gp_values, id, j, detJ,
			dNi_dx, dNi_dy, Me, be);
}

static double get_jacobian(const void *msh, uint32_t id,  double Jinv[4])
{
	/* Compute Jacobian derivatives */
	double dx_dpsi = 0.0;
	double dy_dpsi = 0.0;
	double dx_deta = 0.0;
	double dy_deta = 0.0;

	double dNi_dpsi[3] = {-1, 1, 0};
	double dNi_deta[3] = {-1, 0, 1};
	for (int i = 0; i < 3; i++) {
		uint32_t inode = nb_msh3trg_elem_get_adj(msh, id, i);
		double xi = nb_msh3trg_node_get_x(msh, inode);
		double yi = nb_msh3trg_node_get_y(msh, inode);
		dx_dpsi += dNi_dpsi[i] * xi;
		dx_deta += dNi_deta[i] * xi;
		dy_dpsi += dNi_dpsi[i] * yi;
		dy_deta += dNi_deta[i] * yi;
	}

	/* Compute Jacobian inverse and determinant */
	double detJ = dx_dpsi * dy_deta - dy_dpsi * dx_deta;
      
	Jinv[0] =  dy_deta / detJ;
	Jinv[1] = -dy_dpsi / detJ;
	Jinv[2] = -dx_deta / detJ;
	Jinv[3] =  dx_dpsi / detJ;

	return detJ;
}

static void get_deriv(double Jinv[4], double dNi_dx[3], double dNi_dy[3])
{
	double dNi_dpsi[3] = {-1, 1, 0};
	double dNi_deta[3] = {-1, 0, 1};
	for (uint8_t i = 0; i < 3; i++) {
		dNi_dx[i] = Jinv[0] * dNi_dpsi[i] + Jinv[1] * dNi_deta[i];
		dNi_dy[i] = Jinv[2] * dNi_dpsi[i] + Jinv[3] * dNi_deta[i];
	}
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

void vcn_fem_get_error_on_gpoints(const void *const msh,
				  const vcn_fem_elem_t *const elem,
				  uint32_t N_comp,
				  double* gp_values,
				  double* gp_error)
{
	uint32_t N_nod = nb_msh3trg_get_N_nodes(msh);
	uint32_t N_elem = nb_msh3trg_get_N_elems(msh);

	uint32_t N_gp = vcn_fem_elem_get_N_gpoints(elem);
	memset(gp_error, 0, N_elem * N_gp * sizeof(double));

	/* Interpolate strain to nodes */
	double* nodal_values = calloc(N_comp * N_nod, sizeof(double));
	vcn_fem_interpolate_from_gpoints_to_nodes(msh, elem,
						  N_comp, gp_values,
						  nodal_values);

	for (uint32_t i = 0; i < N_elem; i++) {
		for (uint32_t j = 0; j < N_gp; j++) {
			gp_error[i * N_gp + j] = get_gp_error(i, j, msh,
							      elem, N_comp,
							      gp_values,
							      nodal_values);
		}
	}

	/* Free memory */
	free(nodal_values);
}

static double get_gp_error(uint32_t id_elem, int id_gp,
			   const void *msh,
			   const vcn_fem_elem_t *const elem,
			   uint32_t N_comp,
			   double* gp_values,
			   double* nodal_values)
{
	double *value_gp = alloca(N_comp * sizeof(*value_gp));
	memset(value_gp, 0, N_comp * sizeof(*value_gp));

	uint8_t N_nodes = vcn_fem_elem_get_N_nodes(elem);
	for (uint32_t i = 0; i < N_nodes; i++) {
		uint32_t vid = nb_msh3trg_elem_get_adj(msh, id_elem,
						       N_nodes + i);
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
