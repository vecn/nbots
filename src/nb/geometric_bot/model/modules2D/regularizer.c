#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>

#include "nb/math_bot.h"
#include "nb/eigen_bot.h"
#include "nb/geometric_bot/utils2D.h"

#include "nb/geometric_bot/model/model2D_struct.h"
#include "nb/geometric_bot/model/model2D.h"
#include "nb/geometric_bot/model/modules2D/regularizer.h"

#define GET_1_EDGE_VTX(model, i) ((model)->edge[(i) * 2])
#define GET_2_EDGE_VTX(model, i) ((model)->edge[(i)*2+1])

int vcn_model_regularize(vcn_model_t* model, double lambda,
			 uint32_t N_fixed_vertices,
			 uint32_t* fixed_vertices){
	/* Allocate auxiliar structures */
	double* b = calloc(2 * model->N, sizeof(*b));
    
	/* Allocate sparse matrix */
	vcn_graph_t *graph = vcn_model_get_vtx_graph(model);
	vcn_sparse_t* A = vcn_sparse_create(graph, NULL, 2);
	vcn_graph_destroy(graph);

	/* Assembly system */
	for (uint32_t i = 0; i < model->M; i++) {
		uint32_t id_xi = GET_1_EDGE_VTX(model, i) * 2;
		uint32_t id_yi = GET_1_EDGE_VTX(model, i) * 2+1;
		uint32_t id_xj = GET_2_EDGE_VTX(model, i) * 2;
		uint32_t id_yj = GET_2_EDGE_VTX(model, i) * 2+1;
		/* Assembly matrix 
		 *   _                   _             _  _
		 *  |   1    0   -l   0   |           | xi | 
		 *  |   0    1    0   -l  | = (1 - l) | yi |
		 *  |  -l   0    1    0   |           | xj |
		 *  |_  0   -l   0    1  _|           |_yj_|
		 */

		vcn_sparse_add(A, id_xi, id_xi, 1.0);
		vcn_sparse_add(A, id_xi, id_xj, -lambda);

		vcn_sparse_add(A, id_yi, id_yi, 1.0);
		vcn_sparse_add(A, id_yi, id_yj, -lambda);

		vcn_sparse_add(A, id_xj, id_xi, -lambda);
		vcn_sparse_add(A, id_xj, id_xj, 1.0);

		vcn_sparse_add(A, id_yj, id_yi, -lambda);
		vcn_sparse_add(A, id_yj, id_yj, 1.0);

		/* Assembly RHS vector */
		b[id_xi] += (1.0 - lambda) * model->vertex[id_xi];
		b[id_yi] += (1.0 - lambda) * model->vertex[id_yi];
		b[id_xj] += (1.0 - lambda) * model->vertex[id_xj];
		b[id_yj] += (1.0 - lambda) * model->vertex[id_yj];
	}

	/* Set fixed vertices */
	for (uint32_t i = 0; i < N_fixed_vertices; i++) {
		uint32_t idx = fixed_vertices[i] * 2;
		uint32_t idy = fixed_vertices[i]*2+1;
		vcn_sparse_set_Dirichlet_condition(A, b, idx, model->vertex[idx]);
		vcn_sparse_set_Dirichlet_condition(A, b, idy, model->vertex[idy]);
	}

	/* Solve system */
	uint32_t niter;
	double tol;
	int solver_status = 
		vcn_sparse_solve_CG_precond_Jacobi(A, b, model->vertex, 
						   vcn_sparse_get_size(A),
						   1e-12, &niter, &tol, 1);

	/* Free memory */
	vcn_sparse_destroy(A);
	free(b);

	/* Return status of solver */
	return solver_status;
}
