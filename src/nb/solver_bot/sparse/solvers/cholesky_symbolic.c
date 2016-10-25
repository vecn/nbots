#include <string.h>
#include <stdlib.h>
#include <stdint.h>

#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/solver_bot/sparse/sparse.h"

#include "../sparse_struct.h"

#include "cholesky_symbolic.h"

typedef struct{
	/* Symbolic Cholesky Set
	 * This structure is going to be useful to create
	 * lists that represent sets to implement the 
	 * "symbolic cholesky factorization" algorithm.
	 */
	uint32_t index;
	void *next;
}sc_set;

void nb_sparse_cholesky_symbolic(const nb_sparse_t * const A, 
				 nb_sparse_t *_L, nb_sparse_t* _Lt, 
				 uint32_t ktrunc)
{
	/* The algorithm is going to store in "l" the index of
	 * columns that must have each column of the lower triangular
	 * matrix "L" from A = LL' (without the index of the column "Lii"), 
	 * using only 'ktrunc' diagonals from A.
	 * It could look something like:
	 *   l[0]:   1 -> 2 -> 3 -> 5 -> NULL
	 *   l[1]:   2 -> 4 -> 6 -> 7 -> 8 -> NULL
	 *   l[2]:   3 -> 5 -> 7 -> 9 -> NULL
	 *     :
	 *     :
	 *   l[n-2]: n-1 -> NULL
	 * for a matrix of size n, starting the index from 0.
	 * After each iteration, the algortihm is going to allocate the
	 * i row of the upper matrix "Lt", finally, the matrix "L" will
	 * be allocated in row compress format (nb_sparse_t structure).
	 *
	 * Reference:
	 *  GALLIVAN, Heath, Ng, Ortega, Peyton, Plemmons, Romine, 
	 *    Sameh & Voigt. 
	 *    "Parallel Algorithms for Matrix Computations"
	 *    SIAM 1990 p86-88
	 */
	uint32_t *L_size = nb_allocate_zero_mem(A->N * sizeof(*L_size));

	sc_set** r = nb_allocate_zero_mem(A->N * sizeof(*r));

	for (uint32_t j = 0; j < A->N; j++) { 
		sc_set* lj = NULL;
		uint32_t lj_size = 1;
		uint32_t _i = nb_sparse_bsearch_row(A, j, j, 0,
						 A->rows_size[j]-1);

		/* lj <- aj ************************************************/
		sc_set *iterator_lj;                                     /**/
		for(uint32_t i = _i+1; i<A->rows_size[j]; i++){     /**/
			sc_set* aji = (sc_set*)nb_allocate_mem(sizeof(sc_set));         /**/
			aji->index = A->rows_index[j][i];                      /**/
			aji->next = NULL;                                      /**/
			if(lj == NULL){                                        /**/
				lj = aji;                                            /**/
				iterator_lj = lj;                                    /**/
			}else{                                                 /**/
				iterator_lj->next = aji;                             /**/
				iterator_lj = aji;                                   /**/
			}                                                      /**/
			lj_size++;                                             /**/
		}                                                        /**/
		/***********************************************************/

		sc_set* iterator_rj = r[j];
		while(iterator_rj != NULL){
			/* lj <- lj U {x in li|x-j<k} \ {j} ***********************/
			/*                         ^ This k refers to ktrunc     **/
			/*                           from parameters.            **/
			uint32_t i = iterator_rj->index;                            /**/
			if(lj != NULL){                                         /**/
				for(uint32_t k=1; k<_Lt->rows_size[i]; k++){     /**/
					uint32_t index = _Lt->rows_index[i][k];                 /**/
					if(index != j && index - j < ktrunc){               /**/
						sc_set* next_itr_lj = lj;                         /**/
						sc_set* iterator_lj = NULL;                       /**/
						short flag1 = 0;          /* {index} > all in lj   **/
						while(next_itr_lj != NULL && flag1 == 0){         /**/
							if(next_itr_lj->index > index)                  /**/
								flag1 = 1;            /* {index} < some in lj  **/
							else if(next_itr_lj->index == index)            /**/
								flag1 = 2;            /* {index} in j          **/
							else{                                           /**/
								iterator_lj = next_itr_lj;                    /**/
								next_itr_lj = (sc_set*)next_itr_lj->next;     /**/
							}                                               /**/
						}                                                 /**/
						if(flag1 != 2){                                   /**/
							sc_set* node = (sc_set*)nb_allocate_mem(sizeof(sc_set)); /**/
							node->index = index;                            /**/
							node->next = next_itr_lj;                       /**/
							if(iterator_lj != NULL)                         /**/
								iterator_lj->next = node;                     /**/
							else                                            /**/
								lj = node;                                    /**/
							lj_size++;                                      /**/
						} /* end if (flag1 != 2)                           **/
					}                                                   /**/
				}                                                     /**/
			}else{                                                  /**/
				short flag = 0;             /* {index} == {j}          **/
				uint32_t k = 1;                                           /**/
				while(k < _Lt->rows_size[i] && flag == 0){            /**/
					if(_Lt->rows_index[i][k] != j)                      /**/
						flag = 1;               /* {index} != {j}          **/
					else                                                /**/
						k++;                                              /**/
				}                                                     /**/
				if(flag == 1){                                        /**/
					lj = (sc_set*)nb_allocate_mem(sizeof(sc_set));               /**/
					lj->index = _Lt->rows_index[i][k];                  /**/
					lj->next = NULL;                                    /**/
					lj_size++;                                          /**/
					k++;                                                /**/
				}                                                     /**/
				sc_set* iterator_lj = lj;                             /**/
				while(k < _Lt->rows_size[i]){                         /**/
					if(_Lt->rows_index[i][k] != j){                     /**/
						sc_set* node = (sc_set*)nb_allocate_mem(sizeof(sc_set));   /**/
						node->index = _Lt->rows_index[i][k];              /**/
						node->next = NULL;                                /**/
						iterator_lj->next = node;                         /**/
						lj_size++;                                        /**/
						iterator_lj = (sc_set*)iterator_lj->next;         /**/
					}                                                   /**/
					k++;                                                /**/
				}                                                     /**/
			}                                                       /**/
			/**********************************************************/
			iterator_rj = (sc_set*)iterator_rj->next;
		}

		if(lj != NULL){
			/* p <- min(i in lj) if (lj != NULL) */
			uint32_t p = lj->index;

			/* rp <- rp U {j} *****************************************/
			if(r[p] != NULL){                                       /**/
				sc_set* iterator_rp = r[p];                           /**/
				short flag = 0;                /* {j} is not in rp     **/
				while(iterator_rp->next != NULL && flag == 0){        /**/
					if(iterator_rp->index == j)                         /**/
						flag = 1;                  /* {i} in rp            **/
					iterator_rp = (sc_set*)iterator_rp->next;           /**/
				}                                                     /**/
				if(flag == 0){                                        /**/
					if(iterator_rp->index != j){                        /**/
						sc_set* node = (sc_set*)nb_allocate_mem(sizeof(sc_set));   /**/
						node->index = j;                                  /**/
						node->next = NULL;                                /**/
						iterator_rp->next = node;                         /**/
					}                                                   /**/
				}                                                     /**/
			}else{                                                  /**/
				r[p] = (sc_set*)nb_allocate_mem(sizeof(sc_set));               /**/
				r[p]->index = j;                                      /**/
				r[p]->next = NULL;                                    /**/
			}                                                       /**/
			/**********************************************************/
		}

		/*************** Allocate the jth row of "Lt" **********************/
		_Lt->rows_values[j] = nb_allocate_zero_mem(lj_size *
							   sizeof(double));   /**/
		_Lt->rows_index[j] = nb_allocate_zero_mem(lj_size *
							  sizeof(uint32_t));        /**/
		_Lt->rows_size[j] = lj_size;                                     /**/
		iterator_lj = lj;                                                /**/
		_Lt->rows_index[j][0] = j;                                       /**/
		for(uint32_t k=1; k<lj_size; k++){                          /**/
			_Lt->rows_index[j][k] = iterator_lj->index;                    /**/
			L_size[iterator_lj->index] += 1;                               /**/
			iterator_lj = (sc_set*)iterator_lj->next;                      /**/
		}                                                                /**/
		/* Free memory */                                                /**/
		while(lj != NULL){                                               /**/
			sc_set* lj_free = lj;                                          /**/
			lj = (sc_set*)lj->next;                                        /**/
			nb_free_mem(lj_free);                                                 /**/
		}                                                                /**/
		/*******************************************************************/
    
	}
	/****************************** Allocate "L" ****************************/
	for(uint32_t i=0; i<_L->N; i++){                                 /**/
		L_size[i]++;            /* To include main diagonal */              /**/
		_L->rows_values[i] = nb_allocate_zero_mem(L_size[i] *
							  sizeof(double));     /**/
		_L->rows_index[i] = nb_allocate_zero_mem(L_size[i] *
							 sizeof(uint32_t));          /**/
		_L->rows_size[i] = L_size[i];                                       /**/
	}                                                                     /**/
	/* Cycle to set columns in L */                                       /**/
	uint32_t* L_index = nb_allocate_zero_mem(_L->N * sizeof(uint32_t));                    /**/
	/* This "for" can not be parallelizad to guaranty increasing sort */  /**/
	for(uint32_t i=0; i<_Lt->N; i++){                                /**/
		for(uint32_t j=0; j<_Lt->rows_size[i]; j++){                   /**/
			uint32_t index = _Lt->rows_index[i][j];                               /**/
			_L->rows_index[index][L_index[index]++] = i;                      /**/
		}                                                                   /**/
	}                                                                     /**/
	/************************************************************************/

	/* Free memory */
	nb_free_mem(L_size);
	nb_free_mem(L_index);
	for(uint32_t i=0; i< A->N; i++){
		sc_set* iterator_ri = r[i];
		while(iterator_ri != NULL){
			sc_set* rm = iterator_ri;
			iterator_ri = (sc_set*)iterator_ri->next;
			nb_free_mem(rm);
		}
	}
	nb_free_mem(r);
}
