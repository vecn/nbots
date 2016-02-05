#include <stdio.h>

struct vcn_fem_elem_s {
	uint32_t N_nodes;
	uint32_t N_Gauss_points;
	double *psi, *eta;          /* Normalized space coordinates of Gauss points */
	double * gp_weight;         /* Integration weights of the Gauss points */
	double (**Ni)(double psi, double eta);       /* Shape functions */
	double (**dNi_dpsi)(double psi, double eta); /* Shape functions derivatives */
	double (**dNi_deta)(double psi, double eta); /* Shape functions derivatives */
};

void vcn_fem_elem_destroy(vcn_fem_elem_t* elemtype)
{
	free(elemtype->Ni);
	free(elemtype->dNi_dpsi);
	free(elemtype->dNi_deta);
	free(elemtype->psi);
	free(elemtype->eta);
	free(elemtype->gp_weight);
	free(elemtype);
}

uint32_t vcn_fem_elem_get_N_nodes(const vcn_fem_elem_t *const elemtype)
{
	return elemtype->N_nodes;
}

void* vcn_fem_elem_get_ith_shape_function
		(const vcn_fem_elem_t *const elemtype, uint32_t i)
{
	return elemtype->Ni[i];
}

uint32_t vcn_fem_elem_get_closest_Gauss_Point_to_the_ith_node
		(const vcn_fem_elem_t *const elemtype, uint32_t i)
{
	if (3 == elemtype->N_nodes && 
	    1 == elemtype->N_Gauss_points){
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
