#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "nb/geometric_bot/utils2D.h"
#include "nb/geometric_bot/mesh/elements2D/triangles.h"
#include "nb/pde_bot/finite_element/element.h"
#include "nb/pde_bot/finite_element/gaussp_to_nodes.h"

#include "element_struct.h"

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
	for (uint32_t i = 0; i < N_elements; i++){
		/* Iterate over the nodes of the element */
		for (uint32_t j = 0; j < elemtype->N_nodes; j++){
			uint32_t vj = connectivity_mtx[i*elemtype->N_nodes + j];
			elems_adjacents_to_node[vj*N_room + N_elems_adjacents_to_node[vj]] = i;
			if (N_elems_adjacents_to_node[vj] < N_room - 1)
				N_elems_adjacents_to_node[vj] += 1;
		}
	}
	/* Iterate over vertices to interpolate from the Elemental GP */
	memset(values_interpolated_on_nodes, 0,
	       N_components * N_vertices * sizeof(double));
	for (uint32_t i = 0; i < N_vertices; i++){
		/* Iterate over the elements connected to the vertex */
		double sum_w = 0;
		for (uint32_t k = 0; k < N_elems_adjacents_to_node[i]; k++){
			uint32_t elem_id = elems_adjacents_to_node[i*N_room + k];
			/* Get ID of the vertex relative to the element */
			uint32_t inside_idx = 0;
			while (connectivity_mtx[elem_id * elemtype->N_nodes + inside_idx] != i)
				inside_idx++;

			/* Get id of the closest GP */
			uint32_t id_gp =
				vcn_fem_elem_get_closest_Gauss_Point_to_the_ith_node(elemtype,
										     inside_idx);

			/* Get coordinates to the GP */
			double gp[2] = {0, 0};
			for (uint32_t j = 0; j < elemtype->N_nodes; j++){
				uint32_t vj = connectivity_mtx[k * elemtype->N_nodes + j];
				gp[0] += vertices[vj * 2] *
					elemtype->Ni[j](elemtype->psi[id_gp], elemtype->eta[id_gp]);
				gp[1] += vertices[vj*2+1] *
					elemtype->Ni[j](elemtype->psi[id_gp], elemtype->eta[id_gp]);
			}
			/* Compute weight */
			double wk = 1.0 / vcn_utils2D_get_dist(gp, &(vertices[i*2]));
			sum_w += wk;
			/* Interpolate component by component */
			for (uint32_t c = 0; c < N_components; c++){
				uint32_t cid = elem_id * elemtype->N_Gauss_points + id_gp;
				values_interpolated_on_nodes[i * N_components + c] += 
					wk * values_on_GP_from_elements[cid * N_components + c];
			}
      
		}
		/* Normalize component by component */
		for (uint32_t c = 0; c < N_components; c++)
			values_interpolated_on_nodes[i * N_components + c] /= sum_w;
	}

	/* Free memory */
	free(N_elems_adjacents_to_node);
	free(elems_adjacents_to_node);
}
