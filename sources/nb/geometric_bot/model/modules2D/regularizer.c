#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>

#include "nb/math_bot.h"
#include "nb/memory_bot.h"
#include "nb/solver_bot.h"
#include "nb/geometric_bot/utils2D.h"

#include "nb/geometric_bot/model/model2D_struct.h"
#include "nb/geometric_bot/model/model2D.h"
#include "nb/geometric_bot/model/modules2D/regularizer.h"

#define GET_1_EDGE_VTX(model, i) ((model)->edge[(i) * 2])
#define GET_2_EDGE_VTX(model, i) ((model)->edge[(i)*2+1])

static void init_matrix(const nb_model_t *model, nb_sparse_t **A);
static void assemble_system(const nb_model_t *model, nb_sparse_t *A,
			    double *b, double lambda);

int nb_model_regularize(nb_model_t* model, double lambda,
			 uint32_t N_fixed_vertices,
			 uint32_t* fixed_vertices)
{
	uint32_t memsize = 2 * model->N * sizeof(double);
	char *memblock = nb_soft_allocate_mem(memsize);
	double* b = (void*) memblock;
	memset(b, 0, 2 * model->N  * sizeof(*b));
    
	nb_sparse_t *A;
	init_matrix(model, &A);

	assemble_system(model, A, b, lambda);

	/* Set fixed vertices */
	for (uint32_t i = 0; i < N_fixed_vertices; i++) {
		uint32_t idx = fixed_vertices[i] * 2;
		uint32_t idy = fixed_vertices[i]*2+1;
		nb_sparse_set_Dirichlet_condition(A, b, idx,
						   model->vertex[idx]);
		nb_sparse_set_Dirichlet_condition(A, b, idy,
						   model->vertex[idy]);
	}

	int solver_status = 
		nb_sparse_solve_CG_precond_Jacobi(A, b, model->vertex, 
						   nb_sparse_get_size(A),
						   1e-12, NULL, NULL, 1);

	nb_sparse_destroy(A);
	nb_soft_free_mem(memsize, memblock);

	return solver_status;
}

static void init_matrix(const nb_model_t *model, nb_sparse_t **A)
{
	uint32_t memsize = nb_graph_get_memsize();
	char *memblock = nb_soft_allocate_mem(memsize);
	nb_graph_t *graph = (void*) memblock;
	nb_graph_init(graph);
	
	nb_model_load_vtx_graph(model, graph);
	*A = nb_sparse_create(graph, NULL, 2);

	nb_graph_finish(graph);
	nb_soft_free_mem(memsize, memblock);
}

static void assemble_system(const nb_model_t *model, nb_sparse_t *A,
			    double *b, double lambda)
{
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

		nb_sparse_add(A, id_xi, id_xi, 1.0);
		nb_sparse_add(A, id_xi, id_xj, -lambda);

		nb_sparse_add(A, id_yi, id_yi, 1.0);
		nb_sparse_add(A, id_yi, id_yj, -lambda);

		nb_sparse_add(A, id_xj, id_xi, -lambda);
		nb_sparse_add(A, id_xj, id_xj, 1.0);

		nb_sparse_add(A, id_yj, id_yi, -lambda);
		nb_sparse_add(A, id_yj, id_yj, 1.0);

		/* Assembly RHS vector */
		b[id_xi] += (1.0 - lambda) * model->vertex[id_xi];
		b[id_yi] += (1.0 - lambda) * model->vertex[id_yi];
		b[id_xj] += (1.0 - lambda) * model->vertex[id_xj];
		b[id_yj] += (1.0 - lambda) * model->vertex[id_yj];
	}
}
