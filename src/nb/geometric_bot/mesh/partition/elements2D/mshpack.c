#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <alloca.h>

#include "nb/math_bot.h"
#include "nb/container_bot.h"
#include "nb/eigen_bot.h"
#include "nb/geometric_bot/utils2D.h"
#include "nb/geometric_bot/knn/bins2D.h"
#include "nb/geometric_bot/knn/bins2D_iterator.h"
#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/geometric_bot/mesh/partition/elements2D/mshpack.h"

#include "../../mesh2D_structs.h"

#define POW2(a) ((a)*(a))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

struct nb_mshpack_s {
	uint32_t N_elems;
	double* cen;
	double* radii;
	uint32_t* N_ngb;
	uint32_t** ngb;
};

static void clear_pack(nb_mshpack_t* pack);
static void allocate_mem(nb_mshpack_t *pack, uint32_t N_elems);
static uint32_t mesh_enumerate_input_and_steiner_vtx(vcn_mesh_t *mesh);
static void pack_assemble_adjacencies(const vcn_mesh_t *const mesh,
				       nb_mshpack_t *pack,
				       nb_container_t *segments);
static bool vtx_is_forming_input(const msh_vtx_t *const vtx);
static double pack_optimize_assemble_system(nb_container_t *segments,
					     uint32_t N_spheres,
					     double *Xk,
					     vcn_sparse_t *Hk,
					     double *Bk,
					     double *Xb,
					     double gamma);
static int8_t compare_id(const void* const ptrA,
			 const void* const ptrB);
static void pack_optimize(const vcn_mesh_t *const mesh,
			   nb_mshpack_t *pack, 
			   nb_container_t *segments,
			   double *Xk,
			   double overlapping_factor,
			   uint32_t iterations);
static void pack_update_disks(const vcn_mesh_t *const mesh,
			       nb_mshpack_t *pack,
			       double *Xk);

uint32_t nb_mshpack_get_memsize(void)
{
	return sizeof(nb_mshpack_t);
}

void nb_mshpack_init(void *msh)
{
	memset(msh, 0, nb_mshpack_get_memsize());
}

void nb_mshpack_finish(void *msh)
{
	clear_pack(msh);
}

static void clear_pack(nb_mshpack_t* pack)
{
	if (pack->N_elems > 0) {
		free(pack->cen);
		free(pack->radii);
	}
	if (NULL != pack->ngb) {
		for (uint32_t i=0; i < pack->N_elems; i++)
			free(pack->ngb[i]);
		free(pack->N_ngb);
		free(pack->ngb);
		pack->ngb = NULL;
		pack->N_ngb = NULL;
	}
}

void nb_mshpack_copy(void *msh, const void *mshsrc)
{
	/* PENDING */
}

void nb_mshpack_clear(void *msh)
{
	clear_pack(msh);
	memset(msh, 0, nb_mshpack_get_memsize());
}

uint32_t nb_mshpack_get_N_invtx(const void *msh)
{
	return 0;
}

uint32_t nb_mshpack_get_N_insgm(const void *msh)
{
	return 0;
}

uint32_t nb_mshpack_get_N_nodes(const void *msh)
{
	return 0;
}

uint32_t nb_mshpack_get_N_edges(const void *msh)
{
	return 0;
}

uint32_t nb_mshpack_get_N_elems(const void *msh)
{
	const nb_mshpack_t *pack = msh;
	return pack->N_elems;
}

double nb_mshpack_node_get_x(const void *msh, uint32_t id)
{
	return 0.0;
}

double nb_mshpack_node_get_y(const void *msh, uint32_t id)
{
	return 0.0;
}

uint32_t nb_mshpack_edge_get_1n(const void *msh, uint32_t id)
{
	const nb_mshpack_t *pack = msh;
	return pack->N_elems;
}

uint32_t nb_mshpack_edge_get_2n(const void *msh, uint32_t id)
{
	const nb_mshpack_t *pack = msh;
	return pack->N_elems;
}

double nb_mshpack_elem_get_x(const void *msh, uint32_t id)
{
	const nb_mshpack_t *pack = msh;
	return pack->cen[id * 2];
}

double nb_mshpack_elem_get_y(const void *msh, uint32_t id)
{
	const nb_mshpack_t *pack = msh;
	return pack->cen[id*2+1];
}

double nb_mshpack_elem_get_area(const void *msh, uint32_t id)
{
	double r = nb_mshpack_elem_get_radii(msh, id);
	return NB_PI * POW2(r);
}

double nb_mshpack_elem_face_get_length(const void *msh, 
				       uint32_t elem_id,
				       uint16_t face_id)
{
	uint32_t nj = nb_mshpack_elem_get_ngb(msh, elem_id, face_id);
	double ri = nb_mshpack_elem_get_radii(msh, elem_id);
	double rj = nb_mshpack_elem_get_radii(msh, nj);

	double Ri2 = POW2(ri);
	double Rj2 = POW2(rj);
	double Qij2 = POW2(Ri2 - Rj2);

	const nb_mshpack_t *pack = msh;
	double *s1 = &(pack->cen[elem_id * 2]);
	double *s2 = &(pack->cen[nj * 2]);
	double dij2 = vcn_utils2D_get_dist2(s1, s2);

	return sqrt(3 * Ri2 + Rj2 - dij2 - Qij2 / dij2);
}

double nb_mshpack_elem_get_radii(const void *msh, uint32_t id)
{
	const nb_mshpack_t *pack = msh;
	return pack->radii[id*2+1];
}

uint32_t nb_mshpack_elem_get_N_adj(const void *msh, uint32_t id)
{
	return 0;
}

uint32_t nb_mshpack_elem_get_adj(const void *msh,
				 uint32_t elem_id, uint8_t ngb_id)
{
	const nb_mshpack_t *pack = msh;
	return pack->N_elems;
}

uint32_t nb_mshpack_elem_get_N_ngb(const void *msh, uint32_t id)
{
	const nb_mshpack_t *pack = msh;
	return pack->N_ngb[id];
}

uint32_t nb_mshpack_elem_get_ngb(const void *msh,
				 uint32_t elem_id, uint8_t ngb_id)
{
	const nb_mshpack_t *pack = msh;
	return pack->ngb[elem_id][ngb_id];
}

bool nb_mshpack_elem_has_ngb(const void *msh, uint32_t elem_id,
			     uint16_t ngb_id)
{
	uint32_t N_elems = nb_mshpack_get_N_elems(msh);
	uint32_t id = nb_mshpack_elem_get_ngb(msh, elem_id, ngb_id);
	return id < N_elems;
}

uint32_t nb_mshpack_get_invtx(const void *msh, uint32_t id)
{
	const nb_mshpack_t *pack = msh;
	return pack->N_elems;
}

uint32_t nb_mshpack_insgm_get_N_nodes(const void *msh, uint32_t id)
{
	return 0;
}

uint32_t nb_mshpack_insgm_get_node(const void *msh, uint32_t sgm_id,
				     uint32_t node_id)
{
	const nb_mshpack_t *pack = msh;
	return pack->N_elems;
}

double nb_mshpack_distort_with_field(void *msh,
				     nb_partition_entity field_entity,
				     double *disp,
				     double max_disp)
{
	return 0;/* PENDING */
}

void nb_mshpack_extrapolate_elems_to_nodes(const void *msh, uint8_t N_comp,
					   const double *elem_values,
					   double *nodal_values)
{
	; /* NULL statement */
}

void nb_mshpack_load_elem_graph(const void *msh, nb_graph_t *graph)
{
	;/* NULL statement */
}

void nb_mshpack_load_nodal_graph(const void *msh, nb_graph_t *graph)
{
	;/* NULL statement */
}

void nb_mshpack_load_interelem_graph(const void *msh, nb_graph_t *graph)
{
	;/* PENDING */
}

void nb_mshpack_load_from_mesh_with_overlap(void *msh, nb_mesh_t *mesh,
					    double ov_factor)
{
	nb_mshpack_t *pack = msh;
	uint32_t iterations = 100;
	uint32_t N_elems =
		mesh_enumerate_input_and_steiner_vtx((vcn_mesh_t*)mesh);

	if (N_elems != 0) {
		allocate_mem(pack, N_elems);

		nb_container_t* segments = nb_container_create(NB_QUEUE);

		pack_assemble_adjacencies(mesh, pack, segments);

		double* Xk = malloc(3 * pack->N_elems * sizeof(*Xk));
		pack_optimize(mesh, pack, segments, Xk, ov_factor,
			       iterations);
  
		pack_update_disks(mesh, pack, Xk);
	
		/* Free memory */
		nb_container_set_destroyer(segments, free);
		nb_container_destroy(segments);
		free(Xk);
	}
}

static void allocate_mem(nb_mshpack_t *pack, uint32_t N_elems)
{
	pack->N_elems = N_elems;

	pack->cen = malloc(2 * N_elems * sizeof(*(pack->cen)));
	pack->radii = malloc(N_elems * sizeof(*(pack->radii)));

	pack->N_ngb = calloc(N_elems, sizeof(*(pack->N_ngb)));
	pack->ngb = malloc(N_elems * sizeof(*(pack->ngb)));
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

static void pack_assemble_adjacencies(const vcn_mesh_t *const mesh,
				       nb_mshpack_t *pack,
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
			sgm_struct[2] = pack->N_elems;   /* inner/inner segment */
			
			pack->N_ngb[idx1] += 1;
			pack->N_ngb[idx2] += 1;
		}
	}
	nb_iterator_destroy(iter);

	for (uint32_t i=0; i < pack->N_elems; i++)
		pack->ngb[i] = malloc(pack->N_ngb[i] *
				       sizeof(*(pack->ngb[i])));

	uint32_t* ngb_matrix_next_idx = 
		calloc(pack->N_elems, sizeof(*ngb_matrix_next_idx));

	nb_iterator_t* sgm_iter = nb_iterator_create();
	nb_iterator_set_container(sgm_iter, segments);
	while (nb_iterator_has_more(sgm_iter)) {
		const uint32_t* sgm_struct = nb_iterator_get_next(sgm_iter);    
		if(sgm_struct[2] < pack->N_elems)
			/* Does not consider the nodes on the boundary */
			continue;
    		uint32_t idx1 = sgm_struct[0];
		uint32_t idx2 = sgm_struct[1];
		pack->ngb[idx1][ngb_matrix_next_idx[idx1]] = idx2;
		pack->ngb[idx2][ngb_matrix_next_idx[idx2]] = idx1;
		ngb_matrix_next_idx[idx1] += 1;
		ngb_matrix_next_idx[idx2] += 1;
	}
	nb_iterator_destroy(sgm_iter);
	free(ngb_matrix_next_idx);
}

static bool vtx_is_forming_input(const msh_vtx_t *const vtx)
{
	return mvtx_is_type_origin(vtx, INPUT) ||
		mvtx_is_type_location(vtx, ONSEGMENT);
}

static double pack_optimize_assemble_system(nb_container_t *segments,
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

static int8_t compare_id(const void* const ptrA,
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

static void pack_optimize(const vcn_mesh_t *const mesh,
			   nb_mshpack_t *pack, 
			   nb_container_t *segments,
			   double *Xk,
			   double overlapping_factor,
			   uint32_t iterations)
{
	double gamma = 1 - 0.3 * overlapping_factor;
	/* Allocate Optimization vectors */
	double* hk = calloc(3 * pack->N_elems,  sizeof(*hk));
	double* Bk = malloc(3 * pack->N_elems * sizeof(*Bk));
	/* Allocate boundaries positions */
	uint32_t N_input_vtx = vcn_bins2D_get_length(mesh->ug_vtx) - pack->N_elems;
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

	for (uint32_t i = 0; i < pack->N_elems; i++) {
		/* Initial radii */
		double radii = 0.0;
		for (uint32_t j = 0; j < pack->N_ngb[i]; j++) {
			uint32_t j_id = pack->ngb[i][j];
			radii += (0.5/gamma) * vcn_utils2D_get_dist(&(Xk[i*3]), &(Xk[j_id*3]));
		}
		radii /= pack->N_ngb[i];
		Xk[i*3+2] = radii; /* Initial ratio */
	}

	/* Allocate Hessian as a sparse matrix */
	vcn_graph_t graph;
	graph.N = pack->N_elems;
	graph.N_adj = pack->N_ngb;
	graph.adj = pack->ngb;
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
			memset(Bk, 0, 3 * pack->N_elems * sizeof(double));
			vcn_sparse_make_diagonal(Hk, 1.0);
			global_min =       
				pack_optimize_assemble_system(segments, pack->N_elems,
							       Xk, Hk, Bk, Xb, gamma);
			/* Solve system for qk */
			vcn_sparse_solve_CG_precond_Jacobi(Hk, Bk, hk,
							   vcn_sparse_get_size(Hk) * 2,
							   1e-8, NULL, NULL, 1); 
			/* Xk = Xk + hk */
			for (uint32_t i = 0; i < 3 * pack->N_elems; i++)
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

			printf("min(pack): %e       (%i/%i):(%i/%i)\r", /* TEMPORAL */
			       global_min, optim_k, max_optim_iter,      /* TEMPORAL */
			       super_k, max_super_iter);                 /* TEMPORAL */
			fflush(stdout);                                  /* TEMPORAL */
		}
		printf("                                                 \r"); /* TEMP */
		bool change_adjacencies = false;
		nb_container_t** new_ngb = malloc(pack->N_elems * sizeof(*new_ngb));
		for (uint32_t i = 0; i < pack->N_elems; i++) {
			double gap_factor = 1.5;
			new_ngb[i] = nb_container_create(NB_QUEUE);
			nb_container_set_comparer(new_ngb[i], compare_id);
			/* Add current adjacencies close enough */
			for (uint32_t j = 0; j < pack->N_ngb[i]; j++) {
				uint32_t j_id = pack->ngb[i][j];
				if (vcn_utils2D_get_dist2(&(Xk[i*3]), &(Xk[j_id*3])) <
				    POW2(gap_factor*gamma*(Xk[i*3+2] + Xk[j_id*3+2]))){
					uint32_t* id = malloc(sizeof(*id));
					id[0] = j_id;
					nb_container_insert(new_ngb[i], id);
				} else {
					change_adjacencies = true;
				}
			}
			/* Add adjacencies of adjacent nodes which are close enough */
			for (uint32_t j = 0; j < pack->N_ngb[i]; j++) {
				uint32_t j_id = pack->ngb[i][j];
				for (k = 0; k < pack->N_ngb[j_id]; k++) {
					uint32_t k_id = pack->ngb[j_id][k];
					if(k_id == i)
						continue;
					if (NULL != nb_container_exist(new_ngb[i], &k_id))
						continue;
					if (vcn_utils2D_get_dist2(&(Xk[i*3]), &(Xk[k_id*3])) <
					    POW2(gamma*(Xk[i*3+2] + Xk[k_id*3+2]))) {
						uint32_t* id = malloc(sizeof(*id));
						id[0] = k_id;
						nb_container_insert(new_ngb[i], id);
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
				if (sgm_struct[2] ==  pack->N_elems)
					nb_container_delete(segments, sgm_struct);
			}
			nb_iterator_destroy(sgm_iter);

			/* Reallocate adjacencies and insert new inner/inner segments */
			for (uint32_t i = 0; i < pack->N_elems; i++) {
				free(pack->ngb[i]);
				pack->N_ngb[i] = nb_container_get_length(new_ngb[i]);
				pack->ngb[i] = malloc(pack->N_ngb[i] * sizeof(*(pack->ngb[i])));
				uint32_t j = 0;
				while (nb_container_is_not_empty(new_ngb[i])) {
					uint32_t* id = nb_container_delete_first(new_ngb[i]);
					if (id[0] > i) {
						uint32_t* sgm_struct =
							malloc(3 * sizeof(*sgm_struct));
						sgm_struct[0] = i;
						sgm_struct[1] = id[0];
						sgm_struct[2] = pack->N_elems;
						nb_container_insert(segments, sgm_struct);
					}
					pack->ngb[i][j++] = id[0];
					free(id);
				}
				nb_container_destroy(new_ngb[i]);
			}
			vcn_sparse_destroy(Hk);
			
			vcn_graph_t graph;
			graph.N = pack->N_elems;
			graph.N_adj = pack->N_ngb;
			graph.adj = pack->ngb;
			Hk = vcn_sparse_create(&graph, NULL, 3);
		} else {
			for (uint32_t i = 0; i < pack->N_elems; i++) {
				nb_container_set_destroyer(new_ngb[i], free);
				nb_container_destroy(new_ngb[i]);
			}
		}
		free(new_ngb);
	}

	/* Free memory */
	vcn_sparse_destroy(Hk);
	free(hk);
	free(Bk);
	free(Xb);
}

static void pack_update_disks(const vcn_mesh_t *const mesh,
			       nb_mshpack_t *pack,
			       double *Xk)
{

	vcn_bins2D_iter_t* iter = alloca(vcn_bins2D_iter_get_memsize());
	vcn_bins2D_iter_init(iter);
	vcn_bins2D_iter_set_bins(iter, mesh->ug_vtx);
	while (vcn_bins2D_iter_has_more(iter)) {
		const msh_vtx_t* vtx = vcn_bins2D_iter_get_next(iter);
		if (!vtx_is_forming_input(vtx)) {
			uint32_t id = mvtx_get_id(vtx);
			pack->cen[id * 2] = Xk[id * 3]/mesh->scale + mesh->xdisp;
			pack->cen[id*2+1] = Xk[id*3+1]/mesh->scale + mesh->ydisp;
			pack->radii[id] = Xk[id*3+2] / mesh->scale;
		}
	}
	vcn_bins2D_iter_finish(iter);
}


void nb_mshpack_load_from_mesh(void *msh, nb_mesh_t *mesh)
{
	nb_mshpack_load_from_mesh_with_overlap(msh, mesh, 0.0);
}

void nb_mshpack_get_enveloping_box(const void *msh, double box[4])
{
	const nb_mshpack_t *mshpack = msh;

	box[0] = mshpack->cen[0] - mshpack->radii[0];
	box[1] = mshpack->cen[1] - mshpack->radii[0];
	box[2] = mshpack->cen[0] + mshpack->radii[0];
	box[3] = mshpack->cen[1] + mshpack->radii[0];
	for (uint32_t i = 1; i < mshpack->N_elems; i++) {
		if (box[0] > mshpack->cen[i * 2] - mshpack->radii[i])
			box[0] = mshpack->cen[i * 2] - mshpack->radii[i];
		else if (box[2] < mshpack->cen[i * 2] + mshpack->radii[i]) 
			box[2] = mshpack->cen[i * 2] + mshpack->radii[i];
		if (box[1] > mshpack->cen[i*2+1] - mshpack->radii[i]) 
			box[1] = mshpack->cen[i*2+1] - mshpack->radii[i];
		else if (box[3] < mshpack->cen[i*2+1] + mshpack->radii[i])
			box[3] = mshpack->cen[i*2+1] + mshpack->radii[i];
	}
}

bool nb_mshpack_is_vtx_inside(const void *msh, double x, double y)
{
	return false;/* PENDING */
}

void nb_mshpack_build_model(const void *msh, nb_model_t *model)
{
	;/* PENDING */
}

void nb_mshpack_build_model_disabled_elems(const void *msh,
					   const bool *elems_enabled,
					   nb_model_t *model,
					   uint32_t *N_input_vtx,
					   uint32_t **input_vtx)
{
	;/* PENDING */
}
