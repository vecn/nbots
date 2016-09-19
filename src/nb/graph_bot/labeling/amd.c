/**
 * @brief Labeling used to minimize entries in LU decomposition
 * based on the AMD algorithm.
 * @see P. R. Amestoy, T. A. Davis and I. S. Duff. <b>An Approximate
 * Minimum Degree Ordering Algorithm.</b> Journal of Matrix Analysis
 * and Applications (Vol. 17 1996), SIAM, pages 886-905.
 *
 * @brief Multiple Minimum Degree based on Quotient Graphs. The output
 * labeling reduces the fill-in on LU decompositions.
 * @see A. George and J. W.H. Liu. <b>The evolution of the minimum degree
 * ordering algorithm.</b> SIAM Review (Vol. 31 1989), pages 1-19.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <alloca.h>

#include "nb/container_bot.h"
#include "nb/graph_bot/graph.h"

#define MIN(a,b) (((a)<(b))?(a):(b))

static uint32_t approximate_minimum_degree(int *elem_p, uint32_t label,
					   const nb_graph_t *graph,
					   nb_container_t *minimum_degree,
					   uint32_t* N_adj_elements,
					   void*** adj_elements,
					   nb_container_t** adj_variables,
					   char *flag_elements,
					   uint32_t* N_super_variable,
					   void*** super_variable,
					   uint32_t *perm, uint32_t *iperm);
static void update_vars_i_adj_to_p(int *elem_p,
				   uint32_t* N_adj_elements,
				   void*** adj_elements,
				   nb_container_t** adj_variables);
static void remove_redundant_entries(nb_container_t** adj_variables,
				     int i, int p);
static void process_collisions(nb_container_t *colliding_set,
			       nb_container_t *minimum_degree,
			       uint32_t* N_adj_elements,
			       void*** adj_elements,
			       nb_container_t** adj_variables,
			       uint32_t* N_super_variable,
			       void*** super_variable);
static uint32_t exist(uint32_t N, uint32_t *array, uint32_t val);
static nb_container_t* get_collisions(nb_container_t *hash);
static int8_t compare_degree(const void *const A,
			     const void *const B);

static uint32_t hash_function_supervariables(const void *const);

static uint32_t array_insert(uint32_t N, void*** array, void* val);
static uint32_t array_remove(uint32_t N, void*** array, void* val);

static uint32_t exist(uint32_t N, uint32_t *array, uint32_t val)
{
	/* TEMPORAL Exahustive search */
	for (uint32_t i = 0; i < N; i++) {
		if (array[i] == val)
			return i;
	}
	return N;
}


void nb_graph_labeling_amd(const nb_graph_t *const graph,
			   uint32_t *perm, uint32_t* iperm)
/* Amestoy, Davis & Duff algorithm, 1996 
 *   Variables <- Nodes without label assigned.
 *               (Implemented with a list, faster than the AVL in this case)
 *   Elements  <- Nodes with a label.
 *               (Implemented with an array)
 */
{
	/* Allocate AVL to minimize */
	nb_container_t* minimum_degree = nb_container_create(NB_SORTED);
	nb_container_set_comparer(minimum_degree, compare_degree);
	/* List to store labeled nodes */
	char* flag_elements = calloc(graph->N, sizeof(*flag_elements));

	/* Allocate and initialize nodal data */
	nb_container_t** adj_variables = calloc(graph->N, 
						sizeof(*adj_variables));
	void*** adj_elements = calloc(graph->N, sizeof(*adj_elements));
	uint32_t* N_adj_elements = calloc(graph->N, sizeof(*N_adj_elements));
	int* degree = malloc(2 * graph->N * sizeof(*degree));
	void*** super_variable = malloc(graph->N * sizeof(*super_variable));
	uint32_t* N_super_variable =
		calloc(graph->N, sizeof(*N_super_variable));
	for (uint32_t i = 0; i < graph->N; i++) {
		adj_variables[i] = nb_container_create(NB_QUEUE);
		for (uint32_t j = 0; j < graph->N_adj[i]; j++)
			nb_container_insert(adj_variables[i],
					    &(degree[graph->adj[i][j] * 2]));
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
	memset(perm, 0, graph->N * sizeof(*perm));
	uint32_t label = 0;
	while (nb_container_is_not_empty(minimum_degree)) {
		/* Select variable which minimize the approximated degree */
		int* elem_p = nb_container_delete_first(minimum_degree);
		label = approximate_minimum_degree(elem_p, label,
						   graph, minimum_degree,
						   N_adj_elements,
						   adj_elements,
						   adj_variables,
						   flag_elements,
						   N_super_variable,
						   super_variable,
						   perm, iperm);
	}
	

  	for (uint32_t i = 0; i < graph->N; i++) {
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
}

static uint32_t approximate_minimum_degree(int *elem_p, uint32_t label,
					   const nb_graph_t *graph,
					   nb_container_t *minimum_degree,
					   uint32_t* N_adj_elements,
					   void*** adj_elements,
					   nb_container_t** adj_variables,
					   char *flag_elements,
					   uint32_t* N_super_variable,
					   void*** super_variable,
					   uint32_t *perm, uint32_t *iperm)
{
	/* Add variables from adjacent elements */
	for (uint32_t i = 0; i < N_adj_elements[*elem_p]; i++) {
		int *elem_e = (int*)adj_elements[*elem_p][i];
		nb_iterator_t* subiter = nb_iterator_create();
		nb_iterator_set_container(subiter, adj_variables[*elem_e]);
		while (nb_iterator_has_more(subiter)) {
			int* vtx_i = (int*)nb_iterator_get_next(subiter);
			if (vtx_i == elem_p) 
				continue;
			if (NULL == nb_container_exist(adj_variables[*elem_p], vtx_i))
				nb_container_insert(adj_variables[*elem_p], vtx_i);
		}
		nb_iterator_destroy(subiter);
		nb_container_clear(adj_variables[*elem_e]);
	}

	update_vars_i_adj_to_p(elem_p, N_adj_elements, adj_elements,
			       adj_variables);
    
	/* Compute Le \ Lp for all elements */
	memset(flag_elements, 0, graph->N);
	nb_iterator_t *iter = nb_iterator_create();
	nb_iterator_set_container(iter, adj_variables[*elem_p]);
	while (nb_iterator_has_more(iter)) {
		const int *vtx_i = nb_iterator_get_next(iter);
		for (uint32_t i = 0; i < N_adj_elements[*vtx_i]; i++) {
			int* elem_e = (int*)adj_elements[*vtx_i][i];
			if (!flag_elements[*elem_e]) {
				elem_e[1] = nb_container_get_length(adj_variables[*elem_e]);
				flag_elements[*elem_e] = 1;
			} else {
				elem_e[1] -= N_super_variable[*vtx_i];
			}
		}
	}
	nb_iterator_destroy(iter);

	/* Compute approximate degree */
	iter = nb_iterator_create();
	nb_iterator_set_container(iter, adj_variables[*elem_p]);
	while (nb_iterator_has_more(iter)) {
		int* vtx_i = (int*)nb_iterator_get_next(iter);

		int Ai = nb_container_get_length(adj_variables[*vtx_i]);
		int Lp = nb_container_get_length(adj_variables[*elem_p]) - 1;
			
		uint32_t sum_Le_substract_Lp = 0;
		for (uint32_t j=0; j < N_adj_elements[*vtx_i]; j++) {
			int* elem_e = (int*)adj_elements[*vtx_i][j];
			if (*elem_e != *elem_p)
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
	nb_iterator_set_container(iter, adj_variables[*elem_p]);
	while (nb_iterator_has_more(iter)) {
		const int* vtx_i = nb_iterator_get_next(iter);
		void** var = malloc(4 * sizeof(*var));
		var[0] = (int*) vtx_i;
		var[1] = adj_variables[*vtx_i];
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
		process_collisions(colliding_set, minimum_degree,
				   N_adj_elements, adj_elements,
				   adj_variables,
				   N_super_variable, super_variable);
	}
	nb_iterator_destroy(liter_collision_sets);
		
	nb_container_destroy(hash_collisions);
	nb_container_set_destroyer(hash_supervariables, free);
	nb_container_destroy(hash_supervariables);

	for (uint32_t i = 0; i < N_super_variable[*elem_p]; i++) {
		int *elem_p_i = (int*)super_variable[*elem_p][i];
		perm[*elem_p_i] = label;
		if (iperm != NULL)
			iperm[label] = *elem_p_i;
		label += 1;
		elem_p_i[1] = -1;
	}
	free(adj_elements[*elem_p]);
	return label;
}

static void update_vars_i_adj_to_p(int *elem_p,
				   uint32_t* N_adj_elements,
				   void*** adj_elements,
				   nb_container_t** adj_variables)
{
	nb_iterator_t* iter = alloca(nb_iterator_get_memsize());
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, adj_variables[*elem_p]);
	while (nb_iterator_has_more(iter)) {
		const int *vtx_i = nb_iterator_get_next(iter);
		remove_redundant_entries(adj_variables, *vtx_i, *elem_p);
      
		nb_container_delete(adj_variables[*vtx_i], elem_p);

		/* Element absorption */
		for (uint32_t j = 0; j < N_adj_elements[*elem_p]; j++) {
			int *vtx_j = adj_elements[*elem_p][j];
			N_adj_elements[*vtx_i] =
				array_remove(N_adj_elements[*vtx_i],
					     &(adj_elements[*vtx_i]), vtx_j);
		}

		N_adj_elements[*vtx_i] =
			array_insert(N_adj_elements[*vtx_i], 
				     &(adj_elements[*vtx_i]), elem_p);
	}
	nb_iterator_finish(iter);
}

static void remove_redundant_entries(nb_container_t** adj_variables,
				     int i, int p)
{
	printf("-- SLOW 1\n");/* TEMPORAL AQUI VOY */
	nb_iterator_t* iter = alloca(nb_iterator_get_memsize());
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, adj_variables[p]);
	int TEMP = 0;/* TEMPORAL */
	uint32_t len = nb_container_get_length(adj_variables[p]);/* TEMP */
	while (nb_iterator_has_more(iter)) {
		const int* vtx_j = nb_iterator_get_next(iter);
		nb_container_delete(adj_variables[i], vtx_j);
		TEMP += 1;/* TEMPORAL */
		printf(" %i (%p) ->", *vtx_j, vtx_j);/* TEMPORAL */
		if (TEMP > len) {/* TEMPORAL */
			printf(" Out cycle %i/%i", TEMP, len);/* TEMPORAL */
			break;/* TEMPORAL */
		}/* TEMPORAL */
	}
	nb_iterator_finish(iter);
	printf("\n-- SLOW 2\n");/* TEMPORAL */
	if (TEMP > len)/* TEMPORAL */
		exit(1);/* TEMPORAL */
}

static void process_collisions(nb_container_t *colliding_set,
			       nb_container_t *minimum_degree,
			       uint32_t* N_adj_elements,
			       void*** adj_elements,
			       nb_container_t** adj_variables,
			       uint32_t* N_super_variable,
			       void*** super_variable)
{
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
				if ((nb_container_get_length(adj_variables[*vtx_i]) != 
				     nb_container_get_length(adj_variables[*vtx_j])) ||
				    (N_adj_elements[*vtx_i] != N_adj_elements[*vtx_j])) {
					ij_are_indistinguishable = false;
				} else {
					for (uint32_t i=0; i < N_adj_elements[*vtx_i]; i++) {
						bool is_contained = false;
						for (uint32_t j=0; j < N_adj_elements[*vtx_j]; j++) {
							if (adj_elements[*vtx_i][i] == adj_elements[*vtx_j][j]) {
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
						nb_iterator_set_container(iter_i, adj_variables[*vtx_i]);
						nb_iterator_t* iter_j = nb_iterator_create();
						nb_iterator_set_container(iter_j, adj_variables[*vtx_j]);
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
				for (uint32_t k = 0; k < N_adj_elements[*vtx_j]; k++) {
					int *vtx_k = (int*)adj_elements[*vtx_j][k];
					nb_container_delete(adj_variables[*vtx_k], vtx_j);
					if (NULL == nb_container_exist(adj_variables[*vtx_k], vtx_i))
						nb_container_insert(adj_variables[*vtx_k], vtx_i);
				}
	    
				nb_iterator_t* subiter = 
					nb_iterator_create();
				nb_iterator_set_container(subiter, adj_variables[*vtx_j]);
				while (nb_iterator_has_more(subiter)) {
					int* vtx_k = (int*)nb_iterator_get_next(subiter);
					nb_container_delete(adj_variables[*vtx_k], vtx_j);
					if (NULL == nb_container_exist(adj_variables[*vtx_k], vtx_i))
						nb_container_insert(adj_variables[*vtx_k], vtx_i);
				}
				nb_iterator_destroy(subiter);

				/* Add j to super variable i */
				for (uint32_t i=0; i < N_super_variable[*vtx_j]; i++) {
					int *vtx_k = (int*)super_variable[*vtx_j][i];
					N_super_variable[*vtx_i] = 
						array_insert(N_super_variable[*vtx_i],
							     &(super_variable[*vtx_i]), vtx_k);
				}
						
				nb_container_delete(minimum_degree, vtx_i);
				vtx_i[1] -= N_super_variable[*vtx_j];
				nb_container_insert(minimum_degree, vtx_i);

				nb_container_delete(minimum_degree, vtx_j);

				nb_container_clear(adj_variables[*vtx_j]);
				free(adj_elements[*vtx_j]);

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

void nb_graph_labeling_mmd(const nb_graph_t *const graph,
			   uint32_t *perm, uint32_t* iperm)
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
	memset(perm, 0, graph->N * sizeof(*perm));
	uint32_t label = 0;
	while (nb_container_is_not_empty(minimum_degree)) {
		/* Select variable which minimize the approximated degree */
		int* elem_p = nb_container_delete_first(minimum_degree);

		/* Add variables from adjacent elements */
		for (uint32_t i = 0; i < N_adj_elements[*elem_p]; i++) {
			int *elem_e = (int*)adj_elements[*elem_p][i];
			nb_iterator_t* subiter = nb_iterator_create();
			nb_iterator_set_container(subiter, adj_variables[*elem_e]);
			while (nb_iterator_has_more(subiter)) {
				const int *vtx_i = nb_iterator_get_next(subiter);
				if (vtx_i == elem_p)
					continue;
				if (NULL == nb_container_exist(adj_variables[*elem_p], vtx_i))
					nb_container_insert(adj_variables[*elem_p], vtx_i);
			}
			nb_iterator_destroy(subiter);

			nb_container_clear(adj_variables[*elem_e]);
		}

		/* Update variables i adjacent to p */
		nb_iterator_t* iter = nb_iterator_create();
		nb_iterator_set_container(iter, adj_variables[*elem_p]);
		while (nb_iterator_has_more(iter)) {
			int *vtx_i = (int*) nb_iterator_get_next(iter);
			/* Remove redundant entries */
			nb_iterator_t* subiter = nb_iterator_create();
			nb_iterator_set_container(subiter, adj_variables[*elem_p]);
			while (nb_iterator_has_more(subiter)) {
				int* vtx_j = (int*)nb_iterator_get_next(subiter);
				nb_container_delete(adj_variables[*vtx_i], vtx_j);
			}
			nb_iterator_destroy(subiter);
			nb_container_delete(adj_variables[*vtx_i], elem_p);

			/* Element absorption */
			for (uint32_t j=0; j < N_adj_elements[*elem_p]; j++) {
				int *vtx_j = (int*)adj_elements[*elem_p][j];
				N_adj_elements[*vtx_i] =
					array_remove(N_adj_elements[*vtx_i],
						     &(adj_elements[*vtx_i]), vtx_j);
			}

			N_adj_elements[*vtx_i] =
				array_insert(N_adj_elements[*vtx_i], 
					     &(adj_elements[*vtx_i]), elem_p);

			/* Compute minimum degree */
			int Le = 0;
			if (N_adj_elements[*vtx_i] > 0) {
				nb_container_t* union_of_vars_from_elems_from_i;
				int* first_vtx = (int*)adj_elements[*vtx_i][0];
				union_of_vars_from_elems_from_i = 	  
					nb_container_clone(adj_variables[first_vtx[0]]);
				
				for (uint32_t j=1; j < N_adj_elements[*vtx_i]; j++) {
					int *vtx_j = (int*)adj_elements[*vtx_i][j];
					subiter = nb_iterator_create();
					nb_iterator_set_container(subiter, adj_variables[*vtx_j]);
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
			
			int di = nb_container_get_length(adj_variables[*vtx_i]) + Le;

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
		nb_iterator_set_container(iter, adj_variables[*elem_p]);
		while (nb_iterator_has_more(iter)) {
			const int *vtx_i = nb_iterator_get_next(iter);
			void** var = malloc(4 * sizeof(*var));
			var[0] = (int*) vtx_i;
			var[1] = adj_variables[*vtx_i];
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
						if ((nb_container_get_length(adj_variables[*vtx_i]) != 
						     nb_container_get_length(adj_variables[*vtx_j])) ||
						    (N_adj_elements[*vtx_i] != N_adj_elements[*vtx_j])) {
							ij_are_indistinguishable = false;
						} else {
							for (uint32_t i = 0; i < N_adj_elements[*vtx_i]; i++) {
								bool is_contained = false;
								for (uint32_t j = 0; j < N_adj_elements[*vtx_j]; j++) {
									if (adj_elements[*vtx_i][i] == adj_elements[*vtx_j][j]) {
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
								nb_iterator_set_container(iter_i, adj_variables[*vtx_i]);
								nb_iterator_t* iter_j = nb_iterator_create();
								nb_iterator_set_container(iter_j, adj_variables[*vtx_j]);
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
						for (uint32_t k=0; k < N_adj_elements[*vtx_j]; k++) {
							int *vtx_k = (int*)adj_elements[*vtx_j][k];
							nb_container_delete(adj_variables[*vtx_k], vtx_j);
							if (NULL == nb_container_exist(adj_variables[*vtx_k], vtx_i))
								nb_container_insert(adj_variables[*vtx_k], vtx_i);
						}
						nb_iterator_t* subiter = nb_iterator_create();
						nb_iterator_set_container(subiter, adj_variables[*vtx_j]);
						while (nb_iterator_has_more(subiter)) {
							const int *vtx_k = nb_iterator_get_next(subiter);
							nb_container_delete(adj_variables[*vtx_k], vtx_j);
							if (NULL == nb_container_exist(adj_variables[*vtx_k], vtx_i))
								nb_container_insert(adj_variables[*vtx_k], vtx_i);
						}
						nb_iterator_destroy(subiter);
						/* Add j to super variable i */
						for (uint32_t i = 0; i < N_super_variable[*vtx_j]; i++) {
							int *vtx_k = (int*)super_variable[*vtx_j][i];
							N_super_variable[*vtx_i] = 
								array_insert(N_super_variable[*vtx_i],
									     &(super_variable[*vtx_i]), vtx_k);
						}
						nb_container_delete(minimum_degree, vtx_i);
						vtx_i[1] -= N_super_variable[*vtx_j];
						nb_container_insert(minimum_degree, vtx_i);

						nb_container_delete(minimum_degree, vtx_j);

						nb_container_clear(adj_variables[*vtx_j]);
						free(adj_elements[*vtx_j]);

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

		for (uint32_t i = 0; i < N_super_variable[*elem_p]; i++) {
			int *elem_p_i = (int*)super_variable[*elem_p][i];
			perm[*elem_p_i] = label;
			if (NULL != iperm)
				iperm[label] = *elem_p_i;
			label += 1;
		}
		free(adj_elements[*elem_p]);
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
		hash += *vtx;
	}
	nb_iterator_destroy(iter);
	for (uint32_t k=0; k < N_adj_elements[idx]; k++) {
		int *vtx = (int*)adj_elements[idx][k];
		hash += *vtx;
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
