/******************************************************************************
 *   Graph's Bot: Fast graph theory tools.                                    *
 *   2011-2015 Victor Eduardo Cardoso Nungaray                                *
 *   Twitter: @victore_cardoso                                                *
 *   email: victorc@cimat.mx                                                  *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

#include "nb/container_bot.h"
#include "nb/geometric_bot.h"
#include "nb/eigen_bot.h"
#include "nb/graph_bot/graph.h"

#define MIN(a,b) (((a)<(b))?(a):(b))

static uint32_t get_wij_memsize(nb_graph_t *graph);
static uint32_t exist(uint32_t N, uint32_t *array, uint32_t val);
static void graph_clear(nb_graph_t *graph);
static nb_container_t* get_collisions(nb_container_t *hash);
static int8_t compare_degree(const void *const A,
			     const void *const B);

static int8_t compare_sb(const void *const A, const void *const B,
			 const void *const data);

static uint32_t hash_function_supervariables(const void *const);

static uint32_t array_insert(uint32_t N, void*** array, void* val);
static uint32_t array_remove(uint32_t N, void*** array, void* val);

static vcn_sparse_t* sb_build_laplacian(const vcn_graph_t *const graph);
static uint32_t* sb_partition(const vcn_graph_t *const graph, uint32_t k,
			  uint32_t N_perm, uint32_t *perm);

static uint32_t exist(uint32_t N, uint32_t *array, uint32_t val)
{
	/* TEMPORAL Exahustive search */
	for (uint32_t i = 0; i < N; i++) {
		if (array[i] == val)
			return i;
	}
	return N;
}

uint32_t nb_graph_get_memsize(void)
{
	return sizeof(nb_graph_t);
}

void nb_graph_init(nb_graph_t *graph)
{
	memset(graph, 0, nb_graph_get_memsize());
}

void nb_graph_finish(nb_graph_t *graph)
{
	graph_clear(graph);
}

static void graph_clear(nb_graph_t *graph)
{
	free(graph->N_adj); /* Includes graph->adj in the same memblock */

	if (NULL != graph->wi)
		free(graph->wi);

	if (NULL != graph->wij)
		free(graph->wij);
}

void nb_graph_clear(nb_graph_t *graph)
{
	graph_clear(graph);
	memset(graph, 0, nb_graph_get_memsize());	
}

vcn_graph_t* vcn_graph_create(void)
{
	return calloc(1, sizeof(vcn_graph_t));
}

vcn_graph_t* vcn_graph_get_subgraph(const vcn_graph_t *const graph,
				    uint32_t N_nodes, uint32_t *nodes)
{
	vcn_graph_t *subgraph = calloc(1, sizeof(vcn_graph_t));

	subgraph->N = N_nodes;
	subgraph->N_adj = malloc(N_nodes * sizeof(*(subgraph->N_adj)));
	subgraph->adj = calloc(N_nodes, sizeof(*(subgraph->adj)));

	if (NULL != graph->wi)
		subgraph->wi = malloc(N_nodes * sizeof(*(subgraph->wi)));

	if (NULL != graph->wij)
		subgraph->wij = calloc(N_nodes, sizeof(*(subgraph->wij)));

	for (uint32_t i = 0; i < N_nodes; i++) {
		subgraph->N_adj[i] = 0;
		for (uint32_t j = 0; j < graph->N_adj[nodes[i]]; j++) {
			uint32_t id = exist(N_nodes, nodes, graph->adj[nodes[i]][j]);
			if (id < N_nodes)
				subgraph->N_adj[i] += 1;
		}
		if (subgraph->N_adj[i] > 0) {
			subgraph->adj[i] =
				malloc(subgraph->N_adj[i] * sizeof(*(subgraph->adj[i])));
			if (NULL != graph->wij)
				subgraph->wij[i] = 
					malloc(subgraph->N_adj[i] * sizeof(*(subgraph->wij[i])));
			uint32_t cnt = 0;
			for (uint32_t j = 0; j < graph->N_adj[nodes[i]]; j++) {
				uint32_t id = exist(N_nodes, nodes, graph->adj[nodes[i]][j]);
				if (id < N_nodes) {
					subgraph->adj[i][cnt] = id;
					if (NULL != graph->wij)
						subgraph->wij[i][cnt] = graph->wij[nodes[i]][j];
					cnt += 1;
				}
			}
		}
		if (NULL != graph->wi) 
			subgraph->wi[i] = graph->wi[nodes[i]];
  }

  return subgraph;
}

void vcn_graph_destroy(vcn_graph_t *graph)
{
	nb_graph_finish(graph);
	free(graph);
}

void nb_graph_init_vtx_weights(nb_graph_t *graph)
{
	graph->wi = calloc(graph->N, sizeof(*(graph->wi)));
}

void nb_graph_init_edge_weights(nb_graph_t *graph)
{
	uint32_t memsize = get_wij_memsize(graph);
	char *memblock = calloc(1, memsize);
	graph->wij = (void*) memblock;
	char *block =  memblock + graph->N * sizeof(*(graph->wij));
	for (uint32_t i = 0; i < graph->N; i++) {
		graph->wij[i] = (void*) block;
		block += graph->N_adj[i] * sizeof(**(graph->wij));
	}
}

static uint32_t get_wij_memsize(nb_graph_t *graph)
{
	uint32_t size = graph->N * sizeof(*(graph->wij));
	for (uint32_t i = 0; i < graph->N; i++)
		size += graph->N_adj[i] * sizeof(**(graph->wij));
	return size;
}

void nb_graph_finish_vtx_weights(nb_graph_t *graph)
{
	free(graph->wi);
	graph->wi = NULL;
}

void nb_graph_finish_edge_weights(nb_graph_t *graph)
{
	free(graph->wij);
	graph->wij = NULL;
}

uint32_t* vcn_graph_labeling_amd(const vcn_graph_t *const graph, uint32_t* iperm)
/* Amestoy, Davis & Duff algorithm, 1996 
 *   Variables <- Nodes without label assigned.
 (Implemented with a list, faster than the AVL in this case)
 *   Elements  <- Nodes with a label.
                  (Implemented with an array)
 */
{
	/* Allocate AVL to minimize */
	nb_container_t* minimum_degree = nb_container_create(NB_SORTED);
	nb_container_set_comparer(minimum_degree, compare_degree);
	/* List to store labeled nodes */
	char* flag_elements = calloc(graph->N, sizeof(*flag_elements));

	/* Allocate and initialize nodal data */
	nb_container_t** adj_variables = calloc(graph->N, sizeof(*adj_variables));
	void*** adj_elements = calloc(graph->N, sizeof(*adj_elements));
	uint32_t* N_adj_elements = calloc(graph->N, sizeof(*N_adj_elements));
	int* degree = malloc(2 * graph->N * sizeof(*degree));
	void*** super_variable = malloc(graph->N * sizeof(*super_variable));
	uint32_t* N_super_variable =
		calloc(graph->N, sizeof(*N_super_variable));
	for (uint32_t i = 0; i < graph->N; i++) {
		adj_variables[i] = nb_container_create(NB_QUEUE);
		for (uint32_t j = 0; j < graph->N_adj[i]; j++)
			nb_container_insert(adj_variables[i], &(degree[graph->adj[i][j] * 2]));
		adj_elements[i] = calloc(6, sizeof(*(adj_elements[i])));

		degree[i * 2] = i;
		degree[i*2+1] = nb_container_get_length(adj_variables[i]);

		super_variable[i] = calloc(6, sizeof(*super_variable[i]));
		N_super_variable[i] = 
			array_insert(N_super_variable[i], &(super_variable[i]),
				     &(degree[i*2]));
		nb_container_insert(minimum_degree, &(degree[i*2]));
	}

	/* Approximate Minimum degree */
	uint32_t* perm = calloc(graph->N, sizeof(*perm));
	uint32_t label = 0;
	while (nb_container_is_not_empty(minimum_degree)) {
		/* Select variable which minimize the approximated degree */
		int* elem_p = nb_container_delete_first(minimum_degree);

		/* Add variables from adjacent elements */
		for (uint32_t i=0; i < N_adj_elements[elem_p[0]]; i++) {
			int *elem_e = (int*)adj_elements[elem_p[0]][i];
			nb_iterator_t* subiter = nb_iterator_create();
			nb_iterator_set_container(subiter, adj_variables[elem_e[0]]);
			while (nb_iterator_has_more(subiter)) {
				int* vtx_i = (int*)nb_iterator_get_next(subiter);
				if (vtx_i == elem_p) 
					continue;
				if (NULL == nb_container_exist(adj_variables[elem_p[0]], vtx_i))
					nb_container_insert(adj_variables[elem_p[0]], vtx_i);
			}
			nb_iterator_destroy(subiter);
			nb_container_clear(adj_variables[elem_e[0]]);
		}

		/* Update variables i adjacent to p */
		nb_iterator_t* iter = nb_iterator_create();
		nb_iterator_set_container(iter, adj_variables[elem_p[0]]);
		while (nb_iterator_has_more(iter)) {
			const int *vtx_i = nb_iterator_get_next(iter);
			/* Remove redundant entries */
			nb_iterator_t* subiter = nb_iterator_create();
			nb_iterator_set_container(subiter, adj_variables[elem_p[0]]);
			while (nb_iterator_has_more(subiter)) {
				const int* vtx_j = nb_iterator_get_next(subiter);
				nb_container_delete(adj_variables[vtx_i[0]], vtx_j);
			}
			nb_iterator_destroy(subiter);
      
			nb_container_delete(adj_variables[vtx_i[0]], elem_p);

			/* Element absorption */
			for (uint32_t j = 0; j < N_adj_elements[elem_p[0]]; j++) {
				int *vtx_j = adj_elements[elem_p[0]][j];
				N_adj_elements[vtx_i[0]] =
					array_remove(N_adj_elements[vtx_i[0]],
						     &(adj_elements[vtx_i[0]]), vtx_j);
			}

			N_adj_elements[vtx_i[0]] =
				array_insert(N_adj_elements[vtx_i[0]], 
					     &(adj_elements[vtx_i[0]]), elem_p);
		}
		nb_iterator_destroy(iter);
    
		/* Compute Le \ Lp for all elements */
		memset(flag_elements, 0, graph->N);
		iter = nb_iterator_create();
		nb_iterator_set_container(iter, adj_variables[elem_p[0]]);
		while (nb_iterator_has_more(iter)) {
			const int *vtx_i = nb_iterator_get_next(iter);
			for (uint32_t i = 0; i < N_adj_elements[vtx_i[0]]; i++) {
				int* elem_e = (int*)adj_elements[vtx_i[0]][i];
				if (!flag_elements[elem_e[0]]) {
					elem_e[1] = nb_container_get_length(adj_variables[elem_e[0]]);
					flag_elements[elem_e[0]] = 1;
				} else {
					elem_e[1] -= N_super_variable[vtx_i[0]];
				}
			}
		}
		nb_iterator_destroy(iter);

		/* Compute approximate degree */
		iter = nb_iterator_create();
		nb_iterator_set_container(iter, adj_variables[elem_p[0]]);
		while (nb_iterator_has_more(iter)) {
			int* vtx_i = (int*)nb_iterator_get_next(iter);

			int Ai = nb_container_get_length(adj_variables[vtx_i[0]]);
			int Lp = nb_container_get_length(adj_variables[elem_p[0]]) - 1;
			
			uint32_t sum_Le_substract_Lp = 0;
			for (uint32_t j=0; j < N_adj_elements[vtx_i[0]]; j++) {
				int* elem_e = (int*)adj_elements[vtx_i[0]][j];
				if (elem_e[0] != elem_p[0])
					sum_Le_substract_Lp += elem_e[1];
			}
			int di = MIN(vtx_i[1] + Lp, Ai + Lp + sum_Le_substract_Lp);
			di = MIN(di, nb_container_get_length(minimum_degree));
			/* Update degree */
			if (di != vtx_i[1]) {
				nb_container_delete(minimum_degree, vtx_i);
				vtx_i[1] = di;
				nb_container_insert(minimum_degree, vtx_i);
			}
		}
		nb_iterator_destroy(iter);
		
		/* Supervariables detection */
		nb_container_t* hash_supervariables =
			nb_container_create(NB_HASH);
		nb_container_set_key_generator(hash_supervariables,
						hash_function_supervariables);
		
		iter = nb_iterator_create();
		nb_iterator_set_container(iter, adj_variables[elem_p[0]]);
		while (nb_iterator_has_more(iter)) {
			const int* vtx_i = nb_iterator_get_next(iter);
			void** var = malloc(4 * sizeof(*var));
			var[0] = (int*) vtx_i;
			var[1] = adj_variables[vtx_i[0]];
			var[2] = adj_elements;
			var[3] = N_adj_elements;
			nb_container_insert(hash_supervariables, var);
		}
		nb_iterator_destroy(iter);
		
		nb_container_t* hash_collisions = get_collisions(hash_supervariables);

		nb_iterator_t* liter_collision_sets = 
			nb_iterator_create();
		nb_iterator_set_container(liter_collision_sets, hash_collisions);
		while (nb_iterator_has_more(liter_collision_sets)) {
			nb_container_t* colliding_set = (nb_container_t*)
				nb_iterator_get_next(liter_collision_sets);
			nb_iterator_t* liter_i = nb_iterator_create();
			nb_iterator_set_container(liter_i, colliding_set);
			while (nb_iterator_has_more(liter_i)) {
				void** hash_val_i = (void**) nb_iterator_get_next(liter_i);
				int* vtx_i = (int*) hash_val_i[0];
				nb_iterator_t* liter_j = nb_iterator_clone(liter_i);
				while (nb_iterator_has_more(liter_j)) {
					void** hash_val_j = (void**) nb_iterator_get_next(liter_j);
					int* vtx_j = (int*)hash_val_j[0];
					/* Evaluate if vtx_i and vtx_j are indistinguishable */
					bool ij_are_indistinguishable = true;	
					{
						if ((nb_container_get_length(adj_variables[vtx_i[0]]) != 
						     nb_container_get_length(adj_variables[vtx_j[0]])) ||
						    (N_adj_elements[vtx_i[0]] != N_adj_elements[vtx_j[0]])) {
							ij_are_indistinguishable = false;
						} else {
							for (uint32_t i=0; i < N_adj_elements[vtx_i[0]]; i++) {
								bool is_contained = false;
								for (uint32_t j=0; j < N_adj_elements[vtx_j[0]]; j++) {
									if (adj_elements[vtx_i[0]][i] == adj_elements[vtx_j[0]][j]) {
										is_contained = true;
										break;
									}
								}
								if (!is_contained) {
									ij_are_indistinguishable = false;
									break;
								}
							}
							if (ij_are_indistinguishable) {
								nb_iterator_t* iter_i = nb_iterator_create();
								nb_iterator_set_container(iter_i, adj_variables[vtx_i[0]]);
								nb_iterator_t* iter_j = nb_iterator_create();
								nb_iterator_set_container(iter_j, adj_variables[vtx_j[0]]);
								while (nb_iterator_has_more(iter_i) &&
								       nb_iterator_has_more(iter_j)) {
									const int *ui = nb_iterator_get_next(iter_i);
									const int *uj = nb_iterator_get_next(iter_j);
									if (ui != uj) {
										ij_are_indistinguishable = false;
										break;
									}		
								}
								nb_iterator_destroy(iter_i);
								nb_iterator_destroy(iter_j);
							}
						}
					}/* End of evaluation of indistinguishable variables */

					if (ij_are_indistinguishable) {
						/* Remove references */
						for (uint32_t k = 0; k < N_adj_elements[vtx_j[0]]; k++) {
							int *vtx_k = (int*)adj_elements[vtx_j[0]][k];
							nb_container_delete(adj_variables[vtx_k[0]], vtx_j);
							if (NULL == nb_container_exist(adj_variables[vtx_k[0]], vtx_i))
								nb_container_insert(adj_variables[vtx_k[0]], vtx_i);
						}
	    
						nb_iterator_t* subiter = 
							nb_iterator_create();
						nb_iterator_set_container(subiter, adj_variables[vtx_j[0]]);
						while (nb_iterator_has_more(subiter)) {
							int* vtx_k = (int*)nb_iterator_get_next(subiter);
							nb_container_delete(adj_variables[vtx_k[0]], vtx_j);
							if (NULL == nb_container_exist(adj_variables[vtx_k[0]], vtx_i))
								nb_container_insert(adj_variables[vtx_k[0]], vtx_i);
						}
						nb_iterator_destroy(subiter);

						/* Add j to super variable i */
						for (uint32_t i=0; i < N_super_variable[vtx_j[0]]; i++) {
							int *vtx_k = (int*)super_variable[vtx_j[0]][i];
							N_super_variable[vtx_i[0]] = 
								array_insert(N_super_variable[vtx_i[0]],
									     &(super_variable[vtx_i[0]]), vtx_k);
						}
						
						nb_container_delete(minimum_degree, vtx_i);
						vtx_i[1] -= N_super_variable[vtx_j[0]];
						nb_container_insert(minimum_degree, vtx_i);

						nb_container_delete(minimum_degree, vtx_j);

						nb_container_clear(adj_variables[vtx_j[0]]);
						free(adj_elements[vtx_j[0]]);

						nb_container_delete(colliding_set, hash_val_j);
						nb_iterator_restart(liter_i);
						break;
					}
				}
				nb_iterator_destroy(liter_j);
			}
			nb_iterator_destroy(liter_i);
			nb_container_destroy(colliding_set);
		}
		nb_iterator_destroy(liter_collision_sets);
		
		nb_container_destroy(hash_collisions);
		nb_container_set_destroyer(hash_supervariables, free);
		nb_container_destroy(hash_supervariables);

		for (uint32_t i = 0; i < N_super_variable[elem_p[0]]; i++) {
			int *elem_p_i = (int*)super_variable[elem_p[0]][i];
			perm[elem_p_i[0]] = label;
			if (iperm != NULL)
				iperm[label] = elem_p_i[0];
			label += 1;
			elem_p_i[1] = -1;
		}
    		free(adj_elements[elem_p[0]]);
	}
  	for (uint32_t i=0; i < graph->N; i++) {
		nb_container_destroy(adj_variables[i]);
		free(super_variable[i]);
	}

	/* Free memory */  
	free(adj_variables);
	free(adj_elements);
	free(N_adj_elements);		       
	free(super_variable);
	free(N_super_variable);
	free(degree);
	free(flag_elements);
	nb_container_destroy(minimum_degree);
	return perm;
}

static nb_container_t* get_collisions(nb_container_t *hash)
{
	int8_t status;
	nb_container_t *collisions =
		nb_container_do(hash, "get_collisions", NULL, &status);
	if (0 != status) {
		printf("\nERROR: nb_container_do()\n");
		printf("       using 'get_collisions'\n");
		printf("       exclusive of NB_HASH\n");
		exit(1);
	}
	return collisions;
}

uint32_t* vcn_graph_labeling_mmd(const vcn_graph_t *const graph, uint32_t* iperm)
/* TEMPORAL: Reuse code, very similar to AMD */
{
	/* Allocate AVL to minimize */
	nb_container_t* minimum_degree = nb_container_create(NB_SORTED);
	nb_container_set_comparer(minimum_degree, compare_degree);
	/* List to store labeled nodes */
	char* flag_elements = calloc(graph->N, sizeof(*flag_elements));
  
	/* Allocate and initialize nodal data */
	nb_container_t** adj_variables = calloc(graph->N, sizeof(*adj_variables));
	void*** adj_elements = calloc(graph->N, sizeof(*adj_elements));
	uint32_t* N_adj_elements = calloc(graph->N, sizeof(*N_adj_elements));
	int* degree = malloc(2 * graph->N * sizeof(*degree));
	void*** super_variable = malloc(graph->N * sizeof(*super_variable));
	uint32_t* N_super_variable = calloc(graph->N, sizeof(*N_super_variable));
	for (uint32_t i = 0; i < graph->N; i++) {
		adj_variables[i] = nb_container_create(NB_QUEUE);
		for (uint32_t j = 0; j < graph->N_adj[i]; j++)
			nb_container_insert(adj_variables[i], &(degree[graph->adj[i][j] * 2]));
		adj_elements[i] = calloc(6, sizeof(*adj_elements[i]));

		degree[i * 2] = i;
		degree[i*2+1] = nb_container_get_length(adj_variables[i]);

		super_variable[i] = calloc(6, sizeof(*super_variable[i]));
		N_super_variable[i] = 
			array_insert(N_super_variable[i], &(super_variable[i]), &(degree[i*2]));

		nb_container_insert(minimum_degree, &(degree[i*2]));
	}

	/* Multiple Minimum degree */
	uint32_t* perm = calloc(graph->N, sizeof(*perm));
	uint32_t label = 0;
	while (nb_container_is_not_empty(minimum_degree)) {
		/* Select variable which minimize the approximated degree */
		int* elem_p = nb_container_delete_first(minimum_degree);

		/* Add variables from adjacent elements */
		for (uint32_t i = 0; i < N_adj_elements[elem_p[0]]; i++) {
			int *elem_e = (int*)adj_elements[elem_p[0]][i];
			nb_iterator_t* subiter = nb_iterator_create();
			nb_iterator_set_container(subiter, adj_variables[elem_e[0]]);
			while (nb_iterator_has_more(subiter)) {
				const int *vtx_i = nb_iterator_get_next(subiter);
				if (vtx_i == elem_p)
					continue;
				if (NULL == nb_container_exist(adj_variables[elem_p[0]], vtx_i))
					nb_container_insert(adj_variables[elem_p[0]], vtx_i);
			}
			nb_iterator_destroy(subiter);

			nb_container_clear(adj_variables[elem_e[0]]);
		}

		/* Update variables i adjacent to p */
		nb_iterator_t* iter = nb_iterator_create();
		nb_iterator_set_container(iter, adj_variables[elem_p[0]]);
		while (nb_iterator_has_more(iter)) {
			int *vtx_i = (int*) nb_iterator_get_next(iter);
			/* Remove redundant entries */
			nb_iterator_t* subiter = nb_iterator_create();
			nb_iterator_set_container(subiter, adj_variables[elem_p[0]]);
			while (nb_iterator_has_more(subiter)) {
				int* vtx_j = (int*)nb_iterator_get_next(subiter);
				nb_container_delete(adj_variables[vtx_i[0]], vtx_j);
			}
			nb_iterator_destroy(subiter);
			nb_container_delete(adj_variables[vtx_i[0]], elem_p);

			/* Element absorption */
			for (uint32_t j=0; j < N_adj_elements[elem_p[0]]; j++) {
				int *vtx_j = (int*)adj_elements[elem_p[0]][j];
				N_adj_elements[vtx_i[0]] =
					array_remove(N_adj_elements[vtx_i[0]],
						     &(adj_elements[vtx_i[0]]), vtx_j);
			}

			N_adj_elements[vtx_i[0]] =
				array_insert(N_adj_elements[vtx_i[0]], 
					     &(adj_elements[vtx_i[0]]), elem_p);

			/* Compute minimum degree */
			int Le = 0;
			if (N_adj_elements[vtx_i[0]] > 0) {
				nb_container_t* union_of_vars_from_elems_from_i;
				int* first_vtx = (int*)adj_elements[vtx_i[0]][0];
				union_of_vars_from_elems_from_i = 	  
					nb_container_clone(adj_variables[first_vtx[0]]);
				
				for (uint32_t j=1; j < N_adj_elements[vtx_i[0]]; j++) {
					int *vtx_j = (int*)adj_elements[vtx_i[0]][j];
					subiter = nb_iterator_create();
					nb_iterator_set_container(subiter, adj_variables[vtx_j[0]]);
					while (nb_iterator_has_more(subiter)) {
						int* vtx_k = (int*)nb_iterator_get_next(subiter);
						if (NULL == nb_container_exist(union_of_vars_from_elems_from_i, vtx_k))
							nb_container_insert(union_of_vars_from_elems_from_i, vtx_k);
					}
					nb_iterator_destroy(subiter);
				}
				Le = nb_container_get_length(union_of_vars_from_elems_from_i) - 1;
				nb_container_destroy(union_of_vars_from_elems_from_i);
			}
			
			int di = nb_container_get_length(adj_variables[vtx_i[0]]) + Le;

			/* Update degree */
			if (di != vtx_i[1]) {
				nb_container_delete(minimum_degree, vtx_i);
				vtx_i[1] = di;
				nb_container_insert(minimum_degree, vtx_i);
			}
		}
		nb_iterator_destroy(iter);
		
		/* Supervariables detection */
		nb_container_t* hash_supervariables =
			nb_container_create(NB_HASH);
		nb_container_set_key_generator(hash_supervariables,
						hash_function_supervariables);
		
		iter = nb_iterator_create();
		nb_iterator_set_container(iter, adj_variables[elem_p[0]]);
		while (nb_iterator_has_more(iter)) {
			const int *vtx_i = nb_iterator_get_next(iter);
			void** var = malloc(4 * sizeof(*var));
			var[0] = (int*) vtx_i;
			var[1] = adj_variables[vtx_i[0]];
			var[2] = adj_elements;
			var[3] = N_adj_elements;
			nb_container_insert(hash_supervariables, var);
		}
		nb_iterator_destroy(iter);
		
		nb_container_t* hash_collisions = get_collisions(hash_supervariables);

		nb_iterator_t* liter_collision_sets = 
			nb_iterator_create();
		nb_iterator_set_container(liter_collision_sets, hash_collisions);
		while (nb_iterator_has_more(liter_collision_sets)) {
			nb_container_t* colliding_set = (nb_container_t*)
				nb_iterator_get_next(liter_collision_sets);
			nb_iterator_t* liter_i = nb_iterator_create();
			nb_iterator_set_container(liter_i, colliding_set);
			while (nb_iterator_has_more(liter_i)) {
				void** hash_val_i = (void**)nb_iterator_get_next(liter_i);
				int* vtx_i = (int*)hash_val_i[0];
				nb_iterator_t* liter_j = nb_iterator_clone(liter_i);
				while (nb_iterator_has_more(liter_j)) {
					void** hash_val_j = (void**)nb_iterator_get_next(liter_j);
					int* vtx_j = (int*)hash_val_j[0];
					/* Evaluate if vtx_i and vtx_j are indistinguishable */
					bool ij_are_indistinguishable = true;	
					{
						if ((nb_container_get_length(adj_variables[vtx_i[0]]) != 
						     nb_container_get_length(adj_variables[vtx_j[0]])) ||
						    (N_adj_elements[vtx_i[0]] != N_adj_elements[vtx_j[0]])) {
							ij_are_indistinguishable = false;
						} else {
							for (uint32_t i = 0; i < N_adj_elements[vtx_i[0]]; i++) {
								bool is_contained = false;
								for (uint32_t j = 0; j < N_adj_elements[vtx_j[0]]; j++) {
									if (adj_elements[vtx_i[0]][i] == adj_elements[vtx_j[0]][j]) {
										is_contained = true;
										break;
									}
								}
								if (!is_contained) {
									ij_are_indistinguishable = false;
									break;
								}
							}
							if (ij_are_indistinguishable) {
								nb_iterator_t* iter_i = nb_iterator_create();
								nb_iterator_set_container(iter_i, adj_variables[vtx_i[0]]);
								nb_iterator_t* iter_j = nb_iterator_create();
								nb_iterator_set_container(iter_j, adj_variables[vtx_j[0]]);
								while (nb_iterator_has_more(iter_i) &&
								       nb_iterator_has_more(iter_j)) {
									const int *ui = nb_iterator_get_next(iter_i);
									const int *uj = nb_iterator_get_next(iter_j);
									if (ui != uj) {
										ij_are_indistinguishable = false;
										break;
									}
								}
								nb_iterator_destroy(iter_i);
								nb_iterator_destroy(iter_j);
							}
						}
					}/* End of evaluation of indistinguishable variables */

					if (ij_are_indistinguishable) {
						/* Remove references */
						for (uint32_t k=0; k < N_adj_elements[vtx_j[0]]; k++) {
							int *vtx_k = (int*)adj_elements[vtx_j[0]][k];
							nb_container_delete(adj_variables[vtx_k[0]], vtx_j);
							if (NULL == nb_container_exist(adj_variables[vtx_k[0]], vtx_i))
								nb_container_insert(adj_variables[vtx_k[0]], vtx_i);
						}
						nb_iterator_t* subiter = nb_iterator_create();
						nb_iterator_set_container(subiter, adj_variables[vtx_j[0]]);
						while (nb_iterator_has_more(subiter)) {
							const int *vtx_k = nb_iterator_get_next(subiter);
							nb_container_delete(adj_variables[vtx_k[0]], vtx_j);
							if (NULL == nb_container_exist(adj_variables[vtx_k[0]], vtx_i))
								nb_container_insert(adj_variables[vtx_k[0]], vtx_i);
						}
						nb_iterator_destroy(subiter);
						/* Add j to super variable i */
						for (uint32_t i = 0; i < N_super_variable[vtx_j[0]]; i++) {
							int *vtx_k = (int*)super_variable[vtx_j[0]][i];
							N_super_variable[vtx_i[0]] = 
								array_insert(N_super_variable[vtx_i[0]],
									     &(super_variable[vtx_i[0]]), vtx_k);
						}
						nb_container_delete(minimum_degree, vtx_i);
						vtx_i[1] -= N_super_variable[vtx_j[0]];
						nb_container_insert(minimum_degree, vtx_i);

						nb_container_delete(minimum_degree, vtx_j);

						nb_container_clear(adj_variables[vtx_j[0]]);
						free(adj_elements[vtx_j[0]]);

						nb_container_delete(colliding_set, hash_val_j);
						nb_iterator_restart(liter_i);
						break;
					}
				}
				nb_iterator_destroy(liter_j);
			}
			nb_iterator_destroy(liter_i);
			nb_container_destroy(colliding_set);
		}
		nb_iterator_destroy(liter_collision_sets);
		
		nb_container_destroy(hash_collisions);
		nb_container_set_destroyer(hash_supervariables, free);
		nb_container_destroy(hash_supervariables);

		for (uint32_t i = 0; i < N_super_variable[elem_p[0]]; i++) {
			int *elem_p_i = (int*)super_variable[elem_p[0]][i];
			perm[elem_p_i[0]] = label;
			if (NULL != iperm)
				iperm[label] = elem_p_i[0];
			label += 1;
		}
		free(adj_elements[elem_p[0]]);
	}

	/* Free memory */ 
	for (uint32_t i = 0; i < graph->N; i++) {
		nb_container_destroy(adj_variables[i]);
		free(super_variable[i]);
	}
	free(adj_variables);
	free(adj_elements);
	free(N_adj_elements);	       
	free(super_variable);
	free(N_super_variable);
	free(degree);
	free(flag_elements);
	nb_container_destroy(minimum_degree);
	return perm;
}

static int8_t compare_degree(const void *const A,
			     const void *const B)
{
	const int *const Ai = A;
	const int *const Bi = B;
	if (Ai[1] < Bi[1])
		return -1;
	if (Ai[1] > Bi[1])
		return 1;
	return 0;
}

static uint32_t hash_function_supervariables(const void *const A)
{
	nb_container_t* vars = ((void**)A)[1];
	void*** adj_elements = ((void**)A)[2];
	uint32_t* N_adj_elements = ((void**)A)[3];
	uint32_t idx = ((int*)((void**)A)[0])[0];
	uint32_t hash = 0;
	nb_iterator_t* iter = nb_iterator_create();
	nb_iterator_set_container(iter, vars);
	while (nb_iterator_has_more(iter)) {
		int* vtx = (int*) nb_iterator_get_next(iter);
		hash += vtx[0];
	}
	nb_iterator_destroy(iter);
	for (uint32_t k=0; k < N_adj_elements[idx]; k++) {
		int *vtx = (int*)adj_elements[idx][k];
		hash += vtx[0];
	}
	return hash;
}

static uint32_t array_insert(uint32_t N, void*** array, void* val)
{
	for (uint32_t i = 0; i < N; i++) {
		if (array[0][i] == val)
			return N;
	}
	if (N % 6 == 0 && N > 0) {
		void** new_array = malloc((N+6) * sizeof(*new_array));
		memcpy(new_array, array[0], N * sizeof(*new_array));
		free(array[0]);
		array[0] = new_array;
	}
	array[0][N] = val;
	return N+1;
}

static uint32_t array_remove(uint32_t N, void*** array, void* val)
{
	for (uint32_t i = 0; i < N; i++) {
		if (array[0][i] == val) {
			for (uint32_t j = i; j < N-1; j++)
				array[0][j] = array[0][j+1];
			if ((N-1)%6 == 5) {
				void** new_array = malloc(N * sizeof(*new_array));
				memcpy(new_array, array[0], (N-1) * sizeof(*new_array));
				free(array[0]);
				array[0] = new_array;
			}
			return N-1;
		}
	}
	return N;
}

static vcn_sparse_t* sb_build_laplacian(const vcn_graph_t *const __restrict graph)
{
	vcn_sparse_t* B = vcn_sparse_create(graph, NULL, 1);
	for (uint32_t i = 0; i < graph->N; i++) {
		vcn_sparse_add(B, i, i, 1e-9); /* Damping to the matrix */
		for (uint32_t j = 0; j < graph->N_adj[i]; j++) {
			double wij = 1.0;
			if (graph->wij != NULL)
				wij = graph->wij[i][j];
			vcn_sparse_add(B, i, i, wij);
			vcn_sparse_add(B, i, graph->adj[i][j], -wij);
		}
	}
	return B;
}

static int8_t compare_sb(const void *const A, const void *const B, 
			 const void *const data)
{
	uint32_t id1 = *((uint32_t*)A);
	uint32_t id2 = *((uint32_t*)B);
	double* fiedler = (double*) data;
	if (fiedler[id1] > fiedler[id2])
		return 1;
	if (fiedler[id1] < fiedler[id2]) 
		return -1;
	return 0;
}

static uint32_t* sb_partition(const vcn_graph_t *const __restrict graph, 
			  uint32_t k, uint32_t N_perm, uint32_t *perm)
{
	vcn_sparse_t* B = sb_build_laplacian(graph);

	/* Compute eigenvector with second nondecreasing eigenvalue */
	double* eigen_vals = malloc(2 * sizeof(*eigen_vals));
	double** eigen_vecs = malloc(2 * sizeof(*eigen_vecs));
	eigen_vecs[0] = malloc(graph->N * sizeof(*eigen_vecs[0]));
	eigen_vecs[1] = malloc(graph->N * sizeof(*eigen_vecs[1]));

	int iter; /* TEMPORAL: Unused */
	vcn_sparse_eigen_ipower(B, NB_SOLVER_LUD,
				2, 0.0, eigen_vecs, eigen_vals,
				&iter, 1e-8, 1);
	vcn_sparse_destroy(B);

	uint32_t* p = (uint32_t*) malloc(graph->N * sizeof(uint32_t));
	for (uint32_t i = 0; i < graph->N; i++)
		p[i] = i;

	vcn_qsort_wd(p, graph->N, sizeof(*p), compare_sb, eigen_vecs[0]);

	free(eigen_vals);
	free(eigen_vecs[0]);
	free(eigen_vecs[1]);
	free(eigen_vecs);

	/* Recursion */
	uint32_t *k_points1 = NULL;
	uint32_t *k_points2 = NULL;
	if (k > 3) {
		uint32_t mid = N_perm / 2 + N_perm % 2;

		vcn_graph_t* subgraph = vcn_graph_get_subgraph(graph, mid, p);
		k_points1 = sb_partition(subgraph, k / 2 + k % 2, mid, p);
		free(subgraph);
		
		subgraph = vcn_graph_get_subgraph(graph, N_perm / 2, &(p[mid]));
		k_points2 = sb_partition(subgraph, k / 2, N_perm / 2, &(p[mid]));
		free(subgraph);
	}

	/* Set indices in the correct permutation */
	uint32_t* perm_cpy = malloc(N_perm * sizeof(*perm_cpy));
	memcpy(perm_cpy, perm, N_perm * sizeof(*perm_cpy));
	for (uint32_t i = 0; i < graph->N; i++)
		perm[i] = perm_cpy[p[i]];
	free(p);

	uint32_t* k_points = calloc(k, sizeof(*k_points));
	if (k == 2) {
		k_points[0] = N_perm / 2  + N_perm % 2;
		k_points[1] = graph->N ;
	} else if (k == 3) {
		k_points[0] = N_perm / 3 + (N_perm % 3) / 2;
		k_points[1] = 2 * k_points[0];
		k_points[2] = graph->N ;
	} else {
		for (uint32_t i = 0; i < k / 2 + k % 2; i++)
			k_points[i] = k_points1[i];
		for (uint32_t i = 0; i < k / 2; i++) {
			uint32_t mid = N_perm / 2 + N_perm % 2;
			uint32_t id = i + k / 2 + k % 2;
			k_points[id] = mid + k_points2[i];
		}
		free(k_points1);
		free(k_points2);
	}
	return k_points;
}

uint32_t* vcn_graph_partition_sb(const vcn_graph_t *const graph, uint32_t k)
{
	uint32_t* perm = malloc(graph->N * sizeof(*perm));
	for (uint32_t i = 0; i < graph->N; i++)
		perm[i] = i;

	uint32_t *k_points = sb_partition(graph, k, graph->N, perm);

	uint32_t *part = malloc(graph->N * sizeof(*part));
	uint32_t k_id = 0;
	for (uint32_t i = 0; i < graph->N; i++) {
		if (i >= k_points[k_id])
			k_id += 1;
		part[perm[i]] = k_id;
	}
	free(perm);
	free(k_points);

	return part;
}
