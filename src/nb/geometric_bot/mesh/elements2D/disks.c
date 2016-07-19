#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <alloca.h>

#include "nb/container_bot.h"
#include "nb/eigen_bot.h"
#include "nb/geometric_bot/utils2D.h"
#include "nb/geometric_bot/knn/bins2D.h"
#include "nb/geometric_bot/knn/bins2D_iterator.h"
#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/geometric_bot/mesh/elements2D/disks.h"

#include "../mesh2D_structs.h"

#define POW2(a) ((a)*(a))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

static uint32_t mesh_enumerate_input_and_steiner_vtx(vcn_mesh_t *mesh);

bool vtx_is_forming_input(const msh_vtx_t *const vtx);

static vcn_mshpack_t* spack_create(uint32_t N_spheres)
{
	vcn_mshpack_t* spack = calloc(1, sizeof(*spack));
	spack->N_spheres = N_spheres;

	spack->centers = malloc(2 * N_spheres * sizeof(*(spack->centers)));
	spack->radii = malloc(N_spheres * sizeof(*(spack->radii)));

	spack->N_adj = calloc(N_spheres, sizeof(*(spack->N_adj)));
	spack->adj = malloc(N_spheres * sizeof(*(spack->adj)));

	return spack;
}

static void spack_assemble_adjacencies(const vcn_mesh_t *const mesh,
				       vcn_mshpack_t *spack,
				       nb_container_t *segments)
{
	nb_iterator_t* iter = nb_iterator_create();
	nb_iterator_set_container(iter, mesh->ht_edge);
	while (nb_iterator_has_more(iter)) {
		msh_edge_t* sgm = (msh_edge_t*)nb_iterator_get_next(iter);
		if (vtx_is_forming_input(sgm->v1) != 
		    vtx_is_forming_input(sgm->v2))
			/* A != B : A XOR B (XOR Operator) */
			continue;
		/* Store segments joining inner/inner       segments
		 * and                    boundary/boundary segments.
		 * Discard                inner/boundary    segments.
		 */

		uint32_t idx1 = mvtx_get_id(sgm->v1);
		uint32_t idx2 = mvtx_get_id(sgm->v2);

		if (vtx_is_forming_input(sgm->v1) &&
		    vtx_is_forming_input(sgm->v2)) {
			msh_vtx_t* vtx;
			if (NULL != sgm->t1) {
				vtx = mtrg_get_opposite_vertex(sgm->t1, sgm);
			} else {
				vtx = mtrg_get_opposite_vertex(sgm->t2, sgm);
				idx1 = mvtx_get_id(sgm->v2);
				idx2 = mvtx_get_id(sgm->v1);
			}
			if (vtx_is_forming_input(vtx))
				/* Does not consider the nodes on the boundary */
				continue;

			uint32_t* sgm_struct = malloc(3 * sizeof(*sgm_struct));
			nb_container_insert(segments, sgm_struct);

			uint32_t idx_vtx = mvtx_get_id(vtx);
			sgm_struct[0] = idx1;
			sgm_struct[1] = idx2;
			/* Opposite vtx id (boundary/boundary segment) */
			sgm_struct[2] = idx_vtx;
		} else {
			uint32_t* sgm_struct = malloc(3 * sizeof(*sgm_struct));
			nb_container_insert(segments, sgm_struct);
			
			sgm_struct[0] = idx1;
			sgm_struct[1] = idx2;
			sgm_struct[2] = spack->N_spheres;   /* inner/inner segment */
			
			spack->N_adj[idx1] += 1;
			spack->N_adj[idx2] += 1;
		}
	}
	nb_iterator_destroy(iter);

	for (uint32_t i=0; i < spack->N_spheres; i++)
		spack->adj[i] = malloc(spack->N_adj[i] *
				       sizeof(*(spack->adj[i])));

	uint32_t* adj_matrix_next_idx = 
		calloc(spack->N_spheres, sizeof(*adj_matrix_next_idx));

	nb_iterator_t* sgm_iter = nb_iterator_create();
	nb_iterator_set_container(sgm_iter, segments);
	while (nb_iterator_has_more(sgm_iter)) {
		const uint32_t* sgm_struct = nb_iterator_get_next(sgm_iter);    
		if(sgm_struct[2] < spack->N_spheres)
			/* Does not consider the nodes on the boundary */
			continue;
    		uint32_t idx1 = sgm_struct[0];
		uint32_t idx2 = sgm_struct[1];
		spack->adj[idx1][adj_matrix_next_idx[idx1]] = idx2;
		spack->adj[idx2][adj_matrix_next_idx[idx2]] = idx1;
		adj_matrix_next_idx[idx1] += 1;
		adj_matrix_next_idx[idx2] += 1;
	}
	nb_iterator_destroy(sgm_iter);
	free(adj_matrix_next_idx);
}

bool vtx_is_forming_input(const msh_vtx_t *const vtx)
{
	return mvtx_is_type_origin(vtx, INPUT) ||
		mvtx_is_type_location(vtx, ONSEGMENT);
}

static double spack_optimize_assemble_system(nb_container_t *segments,
					     uint32_t N_spheres,
					     double *Xk,
					     vcn_sparse_t *Hk,
					     double *Bk,
					     double *Xb,
					     double gamma)
{
	double global_min = 0.0;
	double gl[6];    /* Adjacence gradient */
	nb_iterator_t* sgm_iter = nb_iterator_create();
	nb_iterator_set_container(sgm_iter, segments);
	while (nb_iterator_has_more(sgm_iter)) {
		const uint32_t* sgm_struct = nb_iterator_get_next(sgm_iter);
		if (sgm_struct[2] == N_spheres  /* inner/inner segment */) {
			/* Process as inner segment */
			int id1 = sgm_struct[0];
			int id2 = sgm_struct[1];
			/* Compute link function */
			double diff = 
				POW2(Xk[id1 * 3] - Xk[id2 * 3]) +
				POW2(Xk[id1*3+1] - Xk[id2*3+1]) -
				POW2(gamma * (Xk[id1*3+2] + Xk[id2*3+2]));
			double phi_ij =  POW2(diff);
			global_min += phi_ij;
			/* Compute link gradient */
			gl[0] = 4 * diff * (Xk[id1 * 3] - Xk[id2 * 3]);
			gl[1] = 4 * diff * (Xk[id1*3+1] - Xk[id2*3+1]);
			gl[2] = -4 * POW2(gamma) * diff * (Xk[id1*3+2] + Xk[id2*3+2]);
			gl[3] = 4 * diff * (Xk[id2 * 3] - Xk[id1 * 3]);
			gl[4] = 4 * diff * (Xk[id2*3+1] - Xk[id1*3+1]);
			gl[5] = -4 * POW2(gamma) * diff * (Xk[id1*3+2] + Xk[id2*3+2]);
			/* Store global indices */
			int idx[6];
			idx[0] = id1 * 3;
			idx[1] = id1*3+1;
			idx[2] = id1*3+2;
			idx[3] = id2 * 3;
			idx[4] = id2*3+1;
			idx[5] = id2*3+2;
			/* Assembly global gradient */
			for (uint32_t i = 0; i < 6; i++)
				Bk[idx[i]] -= phi_ij * gl[i];
			/* Assembly global Hessian */
			for (uint32_t i = 0; i < 6; i++) {
				for (uint32_t j = 0; j < 6; j++) {
					double val = gl[i]*gl[j];
					vcn_sparse_add(Hk, idx[i], idx[j], val);
				}
			}
		} else if (sgm_struct[2] < N_spheres) {
			/* boundary/boundary segment */
			/* Process as boundary segment */
			double normal[2];
			int k1 = sgm_struct[0];
			int k2 = sgm_struct[1];
			normal[0] = -(Xb[k2*2+1]-Xb[k1*2+1]);
			normal[1] =   Xb[k2 * 2]-Xb[k1 * 2];
			double normalizer = sqrt(POW2(normal[0])+POW2(normal[1]));
			normal[0] /= normalizer;
			normal[1] /= normalizer;
			int id = sgm_struct[2];
			/* Compute boundary function */
			double d_ik = 
				normal[0]*(Xb[k1 * 2] - Xk[id * 3]) + 
				normal[1]*(Xb[k1*2+1] - Xk[id*3+1]);
			double diff = POW2(d_ik) - POW2(Xk[id*3+2]);
			double phi_ik = POW2(diff);
			global_min += phi_ik;
			/* Compute boundary gradient */
			gl[0] = -4 * diff * d_ik * normal[0];
			gl[1] = -4 * diff * d_ik * normal[1];
			gl[2] = -4 * diff * Xk[id*3+2];
			/* Store global indices */
			int idx[3];
			idx[0] = id * 3;
			idx[1] = id*3+1;
			idx[2] = id*3+2;
			/* Assembly global gradient */
			for (uint32_t i = 0; i < 3; i++)
				Bk[idx[i]] -= phi_ik * gl[i];
			/* Assembly global Hessian */
			for (uint32_t i = 0; i < 3; i++) {
				for (uint32_t j = 0; j < 3; j++) {
					double val = gl[i]*gl[j];
					vcn_sparse_add(Hk, idx[i], idx[j], val);
				}
			}
		}
	}
	nb_iterator_destroy(sgm_iter);
	return global_min;
}

static inline int8_t compare_id(const void* const ptrA,
				const void* const ptrB)
{
	uint32_t* A = (uint32_t*) ptrA;
	uint32_t* B = (uint32_t*) ptrB;
	int8_t out;
	if (*A < *B)
		out = -1;
	else if (*A > *B)
		out = 1;
	else
		out = 0;
	return out;
}

static void spack_optimize(const vcn_mesh_t *const mesh,
			   vcn_mshpack_t *spack, 
			   nb_container_t *segments,
			   double *Xk,
			   double overlapping_factor,
			   uint32_t iterations)
{
	double gamma = 1 - 0.3 * overlapping_factor;
	/* Allocate Optimization vectors */
	double* hk = calloc(3 * spack->N_spheres,  sizeof(*hk));
	double* Bk = malloc(3 * spack->N_spheres * sizeof(*Bk));
	/* Allocate boundaries positions */
	uint32_t N_input_vtx = vcn_bins2D_get_length(mesh->ug_vtx) - spack->N_spheres;
	double* Xb =  malloc(N_input_vtx * 2 * sizeof(*Xb));
  
	/***************** Optimize position + radius ***********************/
	/* Initialize solution */
	vcn_bins2D_iter_t* iter = alloca(vcn_bins2D_iter_get_memsize());
	vcn_bins2D_iter_init(iter);
	vcn_bins2D_iter_set_bins(iter, mesh->ug_vtx);
	while (vcn_bins2D_iter_has_more(iter)) {
		const msh_vtx_t* vtx = vcn_bins2D_iter_get_next(iter);
		/* Get ID */
		uint32_t id = mvtx_get_id(vtx);
		if (vtx_is_forming_input(vtx)) {
			Xb[id * 2] = vtx->x[0];
			Xb[id*2+1] = vtx->x[1];
			/* Does not consider nodes in the boundary */
			continue;
		}
		/* Fill 'Xk' */
		Xk[id * 3] = vtx->x[0];
		Xk[id*3+1] = vtx->x[1];
	}
	vcn_bins2D_iter_finish(iter);

	for (uint32_t i = 0; i < spack->N_spheres; i++) {
		/* Initial radii */
		double radii = 0.0;
		for (uint32_t j = 0; j < spack->N_adj[i]; j++) {
			uint32_t j_id = spack->adj[i][j];
			radii += (0.5/gamma) * vcn_utils2D_get_dist(&(Xk[i*3]), &(Xk[j_id*3]));
		}
		radii /= spack->N_adj[i];
		Xk[i*3+2] = radii; /* Initial ratio */
	}

	/* Allocate Hessian as a sparse matrix */
	vcn_graph_t graph;
	graph.N = spack->N_spheres;
	graph.N_adj = spack->N_adj;
	graph.adj = spack->adj;
	vcn_sparse_t* Hk = vcn_sparse_create(&graph, NULL, 3);

	/* Start to minimize contact gaps and intersections 
	   using Newton Method */
	double global_min = 1e100;
	double prev_global_min = 1e100;
	uint32_t super_k = 0;

	/* Start iterations */
	uint32_t max_optim_iter = MIN(iterations, 1);
	/* TEMPORAL: Verify when to update adjacencies (Bean case) */
	uint32_t max_super_iter = MAX(1, iterations/max_optim_iter +
					       ((iterations%max_optim_iter > 0)?1:0));
	while (global_min > NB_GEOMETRIC_TOL && super_k < max_super_iter) {
		super_k ++;
		uint32_t k = 0;
		uint32_t optim_k = 0;
		while (global_min > NB_GEOMETRIC_TOL && optim_k < max_optim_iter) {
			optim_k ++;
			/* Reset global Hessian and Independent vector */
			memset(Bk, 0, 3 * spack->N_spheres * sizeof(double));
			vcn_sparse_make_diagonal(Hk, 1.0);
			global_min =       
				spack_optimize_assemble_system(segments, spack->N_spheres,
							       Xk, Hk, Bk, Xb, gamma);
			/* Solve system for qk */
			vcn_sparse_solve_CG_precond_Jacobi(Hk, Bk, hk,
							   vcn_sparse_get_size(Hk) * 2,
							   1e-8, NULL, NULL, 1); 
			/* Xk = Xk + hk */
			for (uint32_t i = 0; i < 3 * spack->N_spheres; i++)
				/* Compute next step */
				Xk[i] += hk[i];
    			/* Stop criteria */
			if (prev_global_min < global_min)
				k++;
			else if (prev_global_min > global_min){
				k = 0;
				prev_global_min = global_min;
			}
			if (k > 10)
				break;

			printf("min(spack): %e       (%i/%i):(%i/%i)\r", /* TEMPORAL */
			       global_min, optim_k, max_optim_iter,      /* TEMPORAL */
			       super_k, max_super_iter);                 /* TEMPORAL */
			fflush(stdout);                                  /* TEMPORAL */
		}
		printf("                                                 \r"); /* TEMP */
		bool change_adjacencies = false;
		nb_container_t** new_adj = malloc(spack->N_spheres * sizeof(*new_adj));
		for (uint32_t i = 0; i < spack->N_spheres; i++) {
			double gap_factor = 1.5;
			new_adj[i] = nb_container_create(NB_QUEUE);
			nb_container_set_comparer(new_adj[i], compare_id);
			/* Add current adjacencies close enough */
			for (uint32_t j = 0; j < spack->N_adj[i]; j++) {
				uint32_t j_id = spack->adj[i][j];
				if (vcn_utils2D_get_dist2(&(Xk[i*3]), &(Xk[j_id*3])) <
				    POW2(gap_factor*gamma*(Xk[i*3+2] + Xk[j_id*3+2]))){
					uint32_t* id = malloc(sizeof(*id));
					id[0] = j_id;
					nb_container_insert(new_adj[i], id);
				} else {
					change_adjacencies = true;
				}
			}
			/* Add adjacencies of adjacent nodes which are close enough */
			for (uint32_t j = 0; j < spack->N_adj[i]; j++) {
				uint32_t j_id = spack->adj[i][j];
				for (k = 0; k < spack->N_adj[j_id]; k++) {
					uint32_t k_id = spack->adj[j_id][k];
					if(k_id == i)
						continue;
					if (NULL != nb_container_exist(new_adj[i], &k_id))
						continue;
					if (vcn_utils2D_get_dist2(&(Xk[i*3]), &(Xk[k_id*3])) <
					    POW2(gamma*(Xk[i*3+2] + Xk[k_id*3+2]))) {
						uint32_t* id = malloc(sizeof(*id));
						id[0] = k_id;
						nb_container_insert(new_adj[i], id);
						change_adjacencies = true;
					}
				}
			}
		}

		/* Reallocate sparse matrix if a flip has been made */
		if (change_adjacencies) {
			/* Remove inner/inner segments */
			nb_iterator_t* sgm_iter = nb_iterator_create();
			nb_iterator_set_container(sgm_iter, segments);
			while (nb_iterator_has_more(sgm_iter)) {
				uint32_t* sgm_struct = (uint32_t*)nb_iterator_get_next(sgm_iter);
				if (sgm_struct[2] ==  spack->N_spheres)
					nb_container_delete(segments, sgm_struct);
			}
			nb_iterator_destroy(sgm_iter);

			/* Reallocate adjacencies and insert new inner/inner segments */
			for (uint32_t i = 0; i < spack->N_spheres; i++) {
				free(spack->adj[i]);
				spack->N_adj[i] = nb_container_get_length(new_adj[i]);
				spack->adj[i] = malloc(spack->N_adj[i] * sizeof(*(spack->adj[i])));
				uint32_t j = 0;
				while (nb_container_is_not_empty(new_adj[i])) {
					uint32_t* id = nb_container_delete_first(new_adj[i]);
					if (id[0] > i) {
						uint32_t* sgm_struct =
							malloc(3 * sizeof(*sgm_struct));
						sgm_struct[0] = i;
						sgm_struct[1] = id[0];
						sgm_struct[2] = spack->N_spheres;
						nb_container_insert(segments, sgm_struct);
					}
					spack->adj[i][j++] = id[0];
					free(id);
				}
				nb_container_destroy(new_adj[i]);
			}
			vcn_sparse_destroy(Hk);
			
			vcn_graph_t graph;
			graph.N = spack->N_spheres;
			graph.N_adj = spack->N_adj;
			graph.adj = spack->adj;
			Hk = vcn_sparse_create(&graph, NULL, 3);
		} else {
			for (uint32_t i = 0; i < spack->N_spheres; i++) {
				nb_container_set_destroyer(new_adj[i], free);
				nb_container_destroy(new_adj[i]);
			}
		}
		free(new_adj);
	}

	/* Free memory */
	vcn_sparse_destroy(Hk);
	free(hk);
	free(Bk);
	free(Xb);
}

static uint32_t spack_porosity(vcn_mshpack_t *spack,
			   double porosity_factor,
			   bool include_adjacencies)
{
	uint32_t N_porosity_spheres = (uint32_t)
		((1.0 - 0.5 * porosity_factor) * spack->N_spheres);

	if (N_porosity_spheres >= spack->N_spheres)
		return 0;

	uint32_t N_removed_by_porosity = spack->N_spheres - N_porosity_spheres;
	uint32_t id_divisor_porosity = spack->N_spheres / N_removed_by_porosity;

	if (include_adjacencies) {
		uint32_t* N_adj = malloc(N_porosity_spheres * sizeof(*N_adj));
		uint32_t** adj = malloc(N_porosity_spheres * sizeof(*adj));
		for (uint32_t i = 0; i < spack->N_spheres; i++) {
			if (i % id_divisor_porosity != 0 ||
			    i/id_divisor_porosity >= N_removed_by_porosity) {
				uint32_t id = 
					i - MIN(i/id_divisor_porosity+1, N_removed_by_porosity);
				adj[id] = spack->adj[i];
				N_adj[id] = spack->N_adj[i];
				uint32_t j = 0;
				while (j < N_adj[id]) {
					uint32_t j_id = adj[id][j];
					if (j_id % id_divisor_porosity == 0 && 
					    j_id/id_divisor_porosity < N_removed_by_porosity) {
						uint32_t* local_adj = NULL;
						N_adj[id] -= 1;
						if (N_adj[id] > 0) {
							local_adj = malloc(N_adj[id] * sizeof(*local_adj));
							if (j > 0)
								memcpy(local_adj, adj[id], j * sizeof(*local_adj));
							if (j < N_adj[id])
								memcpy(&(local_adj[j]), &(adj[id][j+1]), 
								       (N_adj[id] - j) * sizeof(*local_adj));
						}
						free(adj[id]);
						adj[id] = local_adj;
					} else {
						adj[id][j] = j_id - 
							MIN(j_id/id_divisor_porosity+1, N_removed_by_porosity);
						j++;
					}
				}
			} else {
				free(spack->adj[i]);
			}
		}
		free(spack->adj);
		free(spack->N_adj);
		spack->adj = adj;
		spack->N_adj = N_adj;
	}
	free(spack->centers);
	free(spack->radii);
	spack->N_spheres = N_porosity_spheres;
	spack->centers =
		malloc(2 * N_porosity_spheres * sizeof(*(spack->centers)));
	spack->radii =
		malloc(N_porosity_spheres * sizeof(*spack->radii));
	return N_removed_by_porosity;
}

static void spack_update_disks_porosity(const vcn_mesh_t *const mesh,
					vcn_mshpack_t *spack,
					double *Xk,
					uint32_t N_removed_by_porosity)
{
	uint32_t id_divisor_porosity = spack->N_spheres / N_removed_by_porosity;
	vcn_bins2D_iter_t* iter = alloca(vcn_bins2D_iter_get_memsize());
	vcn_bins2D_iter_init(iter);
	vcn_bins2D_iter_set_bins(iter, mesh->ug_vtx);
	while (vcn_bins2D_iter_has_more(iter)) {
		const msh_vtx_t* vtx = vcn_bins2D_iter_get_next(iter);
		if (!vtx_is_forming_input(vtx)) {
			uint32_t id = mvtx_get_id(vtx);
			if (id % id_divisor_porosity != 0 ||
			    id / id_divisor_porosity >= N_removed_by_porosity) {
				uint32_t id_corrected = id - 
					MIN(id/id_divisor_porosity + 1, N_removed_by_porosity);
				/* Export centroids */
				spack->centers[id_corrected * 2] = 
					Xk[id * 3]/mesh->scale + mesh->xdisp;
				spack->centers[id_corrected*2+1] = 
					Xk[id*3+1]/mesh->scale + mesh->ydisp;
				spack->radii[id_corrected] = Xk[id*3+2] / mesh->scale;
			}
		}
	}
	vcn_bins2D_iter_finish(iter);
}

static void spack_update_disks(const vcn_mesh_t *const mesh,
			       vcn_mshpack_t *spack,
			       double *Xk)
{

	vcn_bins2D_iter_t* iter = alloca(vcn_bins2D_iter_get_memsize());
	vcn_bins2D_iter_init(iter);
	vcn_bins2D_iter_set_bins(iter, mesh->ug_vtx);
	while (vcn_bins2D_iter_has_more(iter)) {
		const msh_vtx_t* vtx = vcn_bins2D_iter_get_next(iter);
		if (!vtx_is_forming_input(vtx)) {
			uint32_t id = mvtx_get_id(vtx);
			spack->centers[id * 2] = Xk[id * 3]/mesh->scale + mesh->xdisp;
			spack->centers[id*2+1] = Xk[id*3+1]/mesh->scale + mesh->ydisp;
			spack->radii[id] = Xk[id*3+2] / mesh->scale;
		}
	}
	vcn_bins2D_iter_finish(iter);
}

vcn_mshpack_t* vcn_mesh_get_mshpack
        (const vcn_mesh_t *const mesh,
	 bool include_adjacencies,
	 uint32_t iterations,
	 double overlapping_factor,  /* Overlapping percentage [0,1] */
	 double porosity_factor,     /* Porosity percentage [0,1] */
	 uint32_t* (*labeling)(const vcn_graph_t *const))
{
	uint32_t N_spheres =  /* Casting mesh to non-const */
		mesh_enumerate_input_and_steiner_vtx((vcn_mesh_t*)mesh);

	if (N_spheres == 0) {
		return NULL;
	}

	vcn_mshpack_t* spack = spack_create(N_spheres);

	nb_container_t* segments = nb_container_create(NB_QUEUE);

	spack_assemble_adjacencies(mesh, spack, segments);

	double* Xk = malloc(3 * spack->N_spheres * sizeof(*Xk));
	spack_optimize(mesh, spack, segments, Xk, overlapping_factor,
		       iterations);
  
	if (!include_adjacencies)
		vcn_mshpack_clear_adj(spack);
  
	uint32_t N_removed_by_porosity =
		spack_porosity(spack, porosity_factor, include_adjacencies);

	if (N_removed_by_porosity > 0)
		spack_update_disks_porosity(mesh, spack, Xk,
					    N_removed_by_porosity);
	else
		spack_update_disks(mesh, spack, Xk);
	
	/* Free memory */
	nb_container_set_destroyer(segments, free);
	nb_container_destroy(segments);
	free(Xk);

	return spack;
}

void vcn_mshpack_clear_adj(vcn_mshpack_t* spack)
{
	if (NULL != spack->adj) {
		for (uint32_t i=0; i < spack->N_spheres; i++)
			free(spack->adj[i]);
		free(spack->N_adj);
		free(spack->adj);
		spack->adj = NULL;
		spack->N_adj = NULL;
	}
}

void vcn_mshpack_destroy(vcn_mshpack_t* spack)
{
	if (spack->N_spheres > 0) {
		free(spack->centers);
		free(spack->radii);
	}
	vcn_mshpack_clear_adj(spack);
	free(spack);
}

static uint32_t mesh_enumerate_input_and_steiner_vtx(vcn_mesh_t *mesh)
{
	uint32_t N_steiner = 0;
	uint32_t N_input = 0;
	vcn_bins2D_iter_t* iter = alloca(vcn_bins2D_iter_get_memsize());
	vcn_bins2D_iter_init(iter);
	vcn_bins2D_iter_set_bins(iter, mesh->ug_vtx);
	while (vcn_bins2D_iter_has_more(iter)) {
		msh_vtx_t* vtx = (msh_vtx_t*) vcn_bins2D_iter_get_next(iter);
		uint32_t id;
		if (!vtx_is_forming_input(vtx))
			id = N_steiner ++; /* Numeration for steiner points */
		else
			id = N_input ++; /* Numeration for fixed nodes in the boundary */
		mvtx_set_id(vtx, id);
	}
	vcn_bins2D_iter_finish(iter);
	return N_steiner;
}
