#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <alloca.h>

#include "nb/memory_bot.h"
#include "nb/geometric_bot/mesh/partition/elements2D/msh3trg.h"

#define INV_3 0.33333333333333333333333333333

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
static void solve_system(const void *msh,
			 uint8_t N_comp, double *M, double *F,
			 double *nodal_values);


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
			    double *M, double *F)
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
	double Jinv[4];
	double detJ = get_jacobian(msh, elem_id, Jinv);

	double dNi_dx[3];
	double dNi_dy[3];
	get_deriv(Jinv, dNi_dx, dNi_dy);

	double wp = 0.5;
	for (uint8_t i = 0; i < 3; i++) {
		uint32_t id = nb_msh3trg_elem_get_adj(msh, elem_id, i);
		double Ni = INV_3;
		for (uint32_t j = 0; j < 3; j++) {
			double Nj = INV_3;
			double integral = Ni * Nj * detJ * wp;
			M[id] += integral;
		}

		for (int c = 0; c < N_comp; c++) {
			double val = elem_values[elem_id * N_comp + c];
			double integral = Ni * val * detJ * wp;
			F[id * N_comp + c] += integral;
		}
	}
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

static void solve_system(const void *msh,
			 uint8_t N_comp, double *M, double *F,
			 double *nodal_values)
{
	uint32_t N_nodes = nb_msh3trg_get_N_nodes(msh);
	for (uint32_t i = 0; i < N_nodes; i++) {
		for (uint8_t j = 0; j < N_comp; j++) {
			uint32_t id = i * N_comp + j;
			nodal_values[id] = F[id] / M[i];
		}
	}
}
