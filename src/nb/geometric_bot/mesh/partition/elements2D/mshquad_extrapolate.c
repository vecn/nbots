#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <alloca.h>

#include "nb/memory_bot.h"
#include "nb/geometric_bot/mesh/partition/elements2D/mshquad.h"

#include "mshquad_struct.h"

#define INV_3 0.33333333333333333333333333333

#define POW2(a) ((a)*(a))

static void assemble_system(const void *msh, 
			    const double* elem_values,
			    uint32_t N_comp,
			    double *M, double *f);

static void trg_assemble_elem(const void *msh,
			      const double *elem_values,
			      uint32_t elem_id,
			      uint8_t N_comp,
			      double *M, double *F);

static double trg_get_jacobian(const void *msh, uint32_t id, double Jinv[4]);
static void trg_get_deriv(double Jinv[4], double dNi_dx[3], double dNi_dy[3]);

static void quad_assemble_elem(const void *msh,
			      const double *elem_values,
			      uint32_t elem_id,
			      uint8_t N_comp,
			      double *M, double *F);

static double quad_get_jacobian(const void *msh, uint32_t id,
				uint8_t gp, double Jinv[4]);
static void quad_get_deriv(uint8_t gp, double Jinv[4], double dNi_dx[4],
			   double dNi_dy[4]);

static void solve_system(const void *msh,
			 uint8_t N_comp, double *M, double *F,
			 double *nodal_values);

static double distort_using_nodal_field(nb_mshquad_t *msh, double *disp,
					double max_disp);
static double get_max_displacement(uint32_t N, double *disp);
static double distort_using_elem_field(nb_mshquad_t *msh, double *disp,
				       double max_disp);


void nb_mshquad_extrapolate_elems_to_nodes(const void *msh, uint8_t N_comp,
					   const double *elem_values,
					   double *nodal_values)
{
	uint32_t N_nodes = nb_mshquad_get_N_nodes(msh);

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
	uint32_t N_elems = nb_mshquad_get_N_elems(msh);
	for (uint32_t i = 0; i < N_elems; i++) {
		if (nb_mshquad_elem_is_quad(msh, i))
			quad_assemble_elem(msh, elem_values,
					   i, N_comp, M, F);
		else
			trg_assemble_elem(msh, elem_values,
					  i, N_comp, M, F);
	}
}

static void trg_assemble_elem(const void *msh,
			      const double *elem_values,
			      uint32_t elem_id,
			      uint8_t N_comp,
			      double *M, double *F)
{
	double Jinv[4];
	double detJ = trg_get_jacobian(msh, elem_id, Jinv);

	double dNi_dx[3];
	double dNi_dy[3];
	trg_get_deriv(Jinv, dNi_dx, dNi_dy);

	double wp = 0.5;
	for (uint8_t i = 0; i < 3; i++) {
		uint32_t id = nb_mshquad_elem_get_adj(msh, elem_id, i);
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

static double trg_get_jacobian(const void *msh, uint32_t id,
			       double Jinv[4])
{
	/* Compute Jacobian derivatives */
	double dx_dpsi = 0.0;
	double dy_dpsi = 0.0;
	double dx_deta = 0.0;
	double dy_deta = 0.0;

	double dNi_dpsi[3] = {-1, 1, 0};
	double dNi_deta[3] = {-1, 0, 1};
	for (int i = 0; i < 3; i++) {
		uint32_t inode = nb_mshquad_elem_get_adj(msh, id, i);
		double xi = nb_mshquad_node_get_x(msh, inode);
		double yi = nb_mshquad_node_get_y(msh, inode);
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

static void trg_get_deriv(double Jinv[4], double dNi_dx[3], double dNi_dy[3])
{
	double dNi_dpsi[3] = {-1, 1, 0};
	double dNi_deta[3] = {-1, 0, 1};
	for (uint8_t i = 0; i < 3; i++) {
		dNi_dx[i] = Jinv[0] * dNi_dpsi[i] + Jinv[1] * dNi_deta[i];
		dNi_dy[i] = Jinv[2] * dNi_dpsi[i] + Jinv[3] * dNi_deta[i];
	}
}

static void quad_assemble_elem(const void *msh,
			       const double *elem_values,
			       uint32_t elem_id, uint8_t N_comp,
			       double *M, double *F)
{
	double wp = 1.0;
	/* Ni(p,n) = (1 + p pi)(1 + n ni)/4 */
	/* GP location = (+-) sqrt(3) / 3 = (+-) 0.577350269190... */
	double Nk[16] = 
		{0.622008467928, 0.166666666667,
		 0.044658198739, 0.166666666667,
		 /*****************************/
		 0.166666666667, 0.622008467928,
		 0.166666666667, 0.044658198739,
		 /*****************************/
		 0.044658198739, 0.166666666667,
		 0.622008467928, 0.166666666667,
		 /*****************************/
		 0.166666666667, 0.044658198739,
		 0.166666666667, 0.622008467928};

	for (uint8_t gp = 0; gp < 4; gp++) {
		double Jinv[4];
		double detJ = quad_get_jacobian(msh, elem_id, gp, Jinv);

		double dNi_dx[4];
		double dNi_dy[4];
		quad_get_deriv(gp, Jinv, dNi_dx, dNi_dy);
	
		for (uint8_t i = 0; i < 4; i++) {
			double Ni = Nk[i * 4 + gp];
			uint32_t id = nb_mshquad_elem_get_adj(msh, elem_id, i);
			for (uint32_t j = 0; j < 4; j++) {
				double Nj = Nk[j * 4 + gp];
				double integral = Ni * Nj * detJ * wp;
				M[id] += integral;
			}

			for (int c = 0; c < N_comp; c++) {
				double w = 0.25; /* Dividing centroid between GP */
				double val = w * elem_values[elem_id * N_comp + c];
				double integral = Ni * val * detJ * wp;
				F[id * N_comp + c] += integral;
			}
		}
	}	
}

static double quad_get_jacobian(const void *msh, uint32_t id,
				uint8_t gp, double Jinv[4])
{
	/* Compute Jacobian derivatives */
	double dx_dpsi = 0.0;
	double dy_dpsi = 0.0;
	double dx_deta = 0.0;
	double dy_deta = 0.0;

	double dNi_dpsi[16] =
		{-0.394337567297, -0.394337567297,
		 -0.105662432703, -0.105662432703,
		 /*******************************/
		 0.394337567297, 0.394337567297,
		 0.105662432703, 0.105662432703,
		 /*******************************/
		 0.105662432703, 0.105662432703,
		 0.394337567297, 0.394337567297,
		 /*******************************/
		 -0.105662432703, -0.105662432703,
		 -0.394337567297, -0.394337567297};
	double dNi_deta[16] =
		{-0.394337567297, -0.105662432703,
		 -0.105662432703, -0.394337567297,
		 /*******************************/
		 -0.105662432703, -0.394337567297,
		 -0.394337567297, -0.105662432703,
		 /*******************************/
		 0.105662432703, 0.394337567297,
		 0.394337567297, 0.105662432703,
		 /*******************************/
		 0.394337567297, 0.105662432703,
		 0.105662432703, 0.394337567297};

	for (int i = 0; i < 4; i++) {
		uint32_t inode = nb_mshquad_elem_get_adj(msh, id, i);
		double xi = nb_mshquad_node_get_x(msh, inode);
		double yi = nb_mshquad_node_get_y(msh, inode);
		for (uint8_t gp = 0; gp < 4; gp++) {
			dx_dpsi += dNi_dpsi[i * 4 + gp] * xi;
			dx_deta += dNi_deta[i * 4 + gp] * xi;
			dy_dpsi += dNi_dpsi[i * 4 + gp] * yi;
			dy_deta += dNi_deta[i * 4 + gp] * yi;
		}
	}

	/* Compute Jacobian inverse and determinant */
	double detJ = dx_dpsi * dy_deta - dy_dpsi * dx_deta;
      
	Jinv[0] =  dy_deta / detJ;
	Jinv[1] = -dy_dpsi / detJ;
	Jinv[2] = -dx_deta / detJ;
	Jinv[3] =  dx_dpsi / detJ;

	return detJ;
}

static void quad_get_deriv(uint8_t gp, double Jinv[4], double dNi_dx[4],
			   double dNi_dy[4])
{
	double dNi_dpsi[16] =
		{-0.394337567297, -0.394337567297,
		 -0.105662432703, -0.105662432703,
		 /*******************************/
		 0.394337567297, 0.394337567297,
		 0.105662432703, 0.105662432703,
		 /*******************************/
		 0.105662432703, 0.105662432703,
		 0.394337567297, 0.394337567297,
		 /*******************************/
		 -0.105662432703, -0.105662432703,
		 -0.394337567297, -0.394337567297};
	double dNi_deta[16] =
		{-0.394337567297, -0.105662432703,
		 -0.105662432703, -0.394337567297,
		 /*******************************/
		 -0.105662432703, -0.394337567297,
		 -0.394337567297, -0.105662432703,
		 /*******************************/
		 0.105662432703, 0.394337567297,
		 0.394337567297, 0.105662432703,
		 /*******************************/
		 0.394337567297, 0.105662432703,
		 0.105662432703, 0.394337567297};
	for (uint8_t i = 0; i < 4; i++) {
		dNi_dx[i] =
			Jinv[0] * dNi_dpsi[i*4+gp] +
			Jinv[1] * dNi_deta[i*4+gp];
		dNi_dy[i] =
			Jinv[2] * dNi_dpsi[i*4+gp] +
			Jinv[3] * dNi_deta[i*4+gp];
	}
}

static void solve_system(const void *msh,
			 uint8_t N_comp, double *M, double *F,
			 double *nodal_values)
{
	uint32_t N_nodes = nb_mshquad_get_N_nodes(msh);
	for (uint32_t i = 0; i < N_nodes; i++) {
		for (uint8_t j = 0; j < N_comp; j++) {
			uint32_t id = i * N_comp + j;
			nodal_values[id] = F[id] / M[i];
		}
	}
}



double nb_mshquad_distort_with_field(void *msh,
				     nb_partition_entity field_entity,
				     double *disp,
				     double max_disp)
{
	double scale = 1.0;
	if (NB_NODE == field_entity)
		scale = distort_using_nodal_field(msh, disp, max_disp);
	else if (NB_ELEMENT == field_entity)
		scale = distort_using_elem_field(msh, disp, max_disp);
	return scale;
}

static double distort_using_nodal_field(nb_mshquad_t *msh, double *disp,
					double max_disp)
{
	uint32_t N_nodes = nb_mshquad_get_N_nodes(msh);
	double max_field_disp = get_max_displacement(N_nodes, disp);
	double scale = max_disp / max_field_disp;
	
	for (uint32_t i = 0; i < 2 * N_nodes; i++)
		msh->nod[i] += disp[i] * scale;

	return scale;
}

static double get_max_displacement(uint32_t N, double *disp)
{
	double max = 0;
	for (uint32_t i = 0; i < N; i++) {
		double disp2 = POW2(disp[i * 2]) + POW2(disp[i*2+1]);
		if (disp2 > max)
			max = disp2;
	}
	return sqrt(max);
}

static double distort_using_elem_field(nb_mshquad_t *msh, double *disp,
				       double max_disp)
{
	uint32_t N = nb_mshquad_get_N_elems(msh);
	uint32_t memsize = 2 * N * sizeof(double);
	double *nodal_disp = NB_SOFT_MALLOC(memsize);
	nb_mshquad_extrapolate_elems_to_nodes(msh, 2, disp, nodal_disp);

	double scale = distort_using_nodal_field(msh, nodal_disp, max_disp);

	NB_SOFT_FREE(memsize, nodal_disp);
	return scale;
}
