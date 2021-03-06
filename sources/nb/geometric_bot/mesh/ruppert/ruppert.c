#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "nb/math_bot.h"
#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/geometric_bot/utils2D.h"
#include "nb/geometric_bot/mesh/tessellator2D.h"

#include "../tessellator2D_structs.h"
#include "ruppert.h"

#define _NB_MAX_LH_TOLERATED (1.0)
#define _NB_MAX_GRADING_RATIO (10.0)
#define _NB_CC_SHELL_UNIT (1e-3)

#define MAX(a,b) (((a)>(b))?(a):(b))
#define NOT_EXIST(container, obj) (!nb_container_exist((container),(obj)))


typedef struct {
	nb_container_t* avl[64];
	uint32_t length;
} hash_trg_t;

static int8_t compare_edge_size(const void *const edge1_ptr,
				const void *const edge2_ptr);
static int8_t compare_trg_attr(const void *const trg1_ptr,
			       const void *const trg2_ptr);
static hash_trg_t* hash_trg_create(void);
static uint32_t hash_trg_length(const hash_trg_t *const htrg);
static void hash_trg_insert(hash_trg_t *const htrg,
			    msh_trg_t *const trg);
static msh_trg_t* hash_trg_remove_first(hash_trg_t *const htrg);
static void hash_trg_remove(hash_trg_t *const htrg,
			    msh_trg_t *const trg);
static void hash_trg_destroy(hash_trg_t *const htrg);

static void reallocate_bins(nb_tessellator2D_t *const restrict mesh);
static bool check_max_vtx(const nb_tessellator2D_t *const mesh);
static bool check_max_trg(const nb_tessellator2D_t *const mesh);
static bool is_encroached(const msh_edge_t *const sgm);
static msh_vtx_t* get_midpoint(nb_tessellator2D_t *mesh,
			       const msh_edge_t *const sgm);
static void concentric_shell(const msh_edge_t *const sgm, msh_vtx_t *v);
static void get_encroached_triangles
                             (const msh_trg_t *const first_trg_to_check,
			      const msh_vtx_t *const v,
			      nb_container_t *encroached_trg);
static void check_encroached_neighbours(const msh_trg_t *trg,
					nb_container_t *encroached_trg,
					nb_container_t *unencroached_trg,
					nb_container_t *processing_trg);
static void remove_encroached_triangles
             (nb_tessellator2D_t *const mesh,
	      msh_trg_t *const first_trg_to_check,
	      const msh_vtx_t *const v,
	      /* v_orfans_reference: NULL to take an arbitrary reference */
	      const msh_vtx_t *const v_orfans_reference,
	      /* big_trg: NULL if not required */
	      nb_container_t *const big_trg,
	      /* poor_quality_trg: NULL if not required */
	      hash_trg_t *const poor_quality_trg,
	      nb_container_t* orfan_vtx);

static msh_vtx_t* retrg_fan_get_next_trg
                          (uint32_t* N_orfan_vtx,
			   msh_vtx_t ** orfan_vtx,
			   const msh_vtx_t *const v1,
			   const msh_vtx_t *const v2);

static void retriangulate_fan
                             (nb_tessellator2D_t *const mesh,
			      const msh_vtx_t *const v_pivot,
			      const msh_vtx_t *const v_start,
			      nb_container_t *orfan_vtx,
			      /* l_new_trg: NULL if not required */
			      nb_container_t *const l_new_trg);
static void insert_vertex(nb_tessellator2D_t *mesh,
			  msh_trg_t *const trg_containing_cc,
			  const msh_vtx_t *const cc,
			  /* big_trg: NULL if not required */
			  nb_container_t *const big_trg,
			  /* poor_quality_trg: NULL if not required */
			  hash_trg_t *const poor_quality_trg,
			  /* l_new_trg: NULL if not required */
			  nb_container_t *const l_new_trg);

static void insert_midpoint(nb_tessellator2D_t *const mesh,
			    msh_edge_t *const sgm,
			    const msh_vtx_t *const v,
			    /* big_trg: NULL if not required */
			    nb_container_t *const big_trg,
			    /* poor_quality_trg: NULL if not required */
			    hash_trg_t *const poor_quality_trg,
			    msh_edge_t* subsgm[2],
			    /* l_new_trg: NULL if not required */
			    nb_container_t *const l_new_trg);

static void split_encroached_segments(nb_tessellator2D_t *const mesh,
				      nb_container_t *const encroached_sgm,
				      nb_container_t *const big_trg,
				      hash_trg_t *const poor_quality_trg);
						   
static void verify_new_encroachments
                             (nb_tessellator2D_t *const mesh,
			      const msh_vtx_t *const  v,
			      nb_container_t *const l_new_trg,
			      nb_container_t *const encroached_sgm,
			      nb_container_t *const big_trg,
			      hash_trg_t *const poor_quality_trg);

static void initialize_encroached_sgm
                          (nb_tessellator2D_t *const mesh,
			   nb_container_t *const encroached_sgm);

static void check_trg(msh_trg_t *const trg, nb_tessellator2D_t *const mesh,
			   nb_container_t *const big_trg,
			   hash_trg_t *const poor_quality_trg);

static void initialize_big_and_poor_quality_trg
                          (nb_tessellator2D_t *const mesh,
			   nb_container_t *const big_trg,
			   hash_trg_t *const poor_quality_trg);

static msh_trg_t* get_trg_containing_circumcenter
                          (const nb_tessellator2D_t *const mesh,
			   const msh_trg_t *const trg,
			   const msh_vtx_t *const cc);

static void get_sgm_encroached_by_vertex
                          (const nb_tessellator2D_t *const mesh,
			   const msh_trg_t *const trg_containing_vtx,
			   const msh_vtx_t *const vtx,
			   nb_container_t *encroached_sgm);

static bool split_is_permitted(const msh_edge_t *const sgm,
			       double d);

static void get_subsgm_cluster(const msh_edge_t *const sgm,
			       const msh_vtx_t *const vtx,
			       double* smallest_angle,
			       nb_container_t* cluster);

static void delete_bad_trg(nb_tessellator2D_t *mesh,
			   nb_container_t *encroached_sgm,
			   nb_container_t *big_trg,
			   hash_trg_t *poor_quality_trg);
static void delete_poor_quality_trg(nb_tessellator2D_t *mesh,
				    nb_container_t *encroached_sgm,
				    nb_container_t *big_trg,
				    hash_trg_t *poor_quality_trg,
				    msh_trg_t *trg);
static void circumcenter_accepted(nb_tessellator2D_t *mesh,
				  nb_container_t *big_trg,
				  hash_trg_t *poor_quality_trg,
				  msh_trg_t *trg_containing_cc,
				  msh_vtx_t *cc);
static void circumcenter_rejected(nb_tessellator2D_t *mesh,
				  nb_container_t *encroached_sgm,
				  nb_container_t *big_trg,
				  hash_trg_t *poor_quality_trg,
				  nb_container_t *sgm_encroached_by_cc,
				  msh_vtx_t *cc, msh_trg_t *trg);
static bool edge_violates_constrain(const nb_tessellator2D_t *const restrict mesh,
				    const msh_edge_t *const restrict sgm,
				    double *big_ratio);
static bool edge_greater_than_density(const nb_tessellator2D_t *const restrict mesh,
				      const msh_edge_t *const restrict sgm,
				      double *big_ratio);
static double calculate_lh(const nb_tessellator2D_t *const restrict mesh,
			   const msh_edge_t *const restrict sgm,
			   double dist, int N_trapezoids);
static double calculate_h(const nb_tessellator2D_t *const restrict mesh,
			  const msh_edge_t *const restrict sgm,
			  double t);
static bool mtrg_is_too_big(const nb_tessellator2D_t *const restrict mesh,
			    const msh_trg_t *const restrict trg,
			    /* big_ratio could be NULL if not required */
			    double *big_ratio);

void nb_ruppert_refine(nb_tessellator2D_t *restrict mesh)
{
	/* Allocate data structures to allocate encroached elements */
	nb_container_t *encroached_sgm =
		nb_allocate_on_stack(nb_container_get_memsize(NB_SORTED));
	nb_container_init(encroached_sgm, NB_SORTED);
	nb_container_set_comparer(encroached_sgm, compare_edge_size);

	nb_container_t *big_trg =
		nb_allocate_on_stack(nb_container_get_memsize(NB_SORTED));
	nb_container_init(big_trg, NB_SORTED);
	nb_container_set_comparer(big_trg, compare_trg_attr);

	hash_trg_t *poor_quality_trg = hash_trg_create();

	/* Initialize container with encroached segments */
	initialize_encroached_sgm(mesh, encroached_sgm);

	/* Calculate max circumradius to shortest edge ratio allowed */
	if (mesh_get_min_angle(mesh) > NB_MESH_MAX_ANGLE) {
		if (0 == mesh->max_vtx) {
			printf("WARNING in nb_tessellator2D_refine(): ");
			printf("Setting max_vtx = 1000000 to");
			printf("warranty finish.\n");
			/* Set a max because there aren't guarantees to finish */
			mesh->max_vtx = 1000000;
		}
	}

	/* Split encroached segments */
	split_encroached_segments(mesh, encroached_sgm,
				  big_trg, poor_quality_trg);

	/* Initialize hash table with the poor quality triangles */
	initialize_big_and_poor_quality_trg(mesh, big_trg,
					    poor_quality_trg);
	delete_bad_trg(mesh, encroached_sgm, big_trg, poor_quality_trg);

	/* Free data structures */
	nb_container_finish(encroached_sgm);
	nb_container_finish(big_trg);
	hash_trg_destroy(poor_quality_trg);
}

static int8_t compare_edge_size(const void *const edge1_ptr,
				const void *const edge2_ptr)
{
	const msh_edge_t *const edge1 = edge1_ptr;
	const msh_edge_t *const edge2 = edge2_ptr;
	uint32_t a1 = (uint32_t)
		(nb_utils2D_get_dist2(edge1->v1->x, edge1->v2->x) * 1e6);
	uint32_t a2 = (uint32_t)
		(nb_utils2D_get_dist2(edge2->v1->x, edge2->v2->x) * 1e6);
	int8_t out;
	if (a1 > a2) {
		out = -1;
	} else if (a1 < a2) {
		out = 1;
	} else {
		if (edge1 < edge2)
			out = -1;
		else if (edge1 > edge2)
			out = 1;
		else
			out = 0;
	}
	return out;
}

static int8_t compare_trg_attr(const void *const trg1_ptr,
			       const void *const trg2_ptr)
{
	const msh_trg_t *const trg1 = trg1_ptr;
	const msh_trg_t *const trg2 = trg2_ptr;
	uint32_t a1 = (uint32_t)(1e6 * trg1->feature);
	uint32_t a2 = (uint32_t)(1e6 * trg2->feature);
	int8_t out;
	if (a1 > a2) {
		out = -1;
	} else if (a1 < a2) {
		out = 1;
	} else {
		if (trg1 < trg2)
			out = -1;
		else if (trg1 > trg2)
			out = 1;
		else
			out = 0;
	}
	return out;
}

bool nb_ruppert_insert_vtx(nb_tessellator2D_t *mesh, const double vertex[2])
{
	msh_vtx_t* new_vtx = mvtx_create(mesh);

	/* Move and scale new vertex */
	new_vtx->x[0] = mesh->scale * (vertex[0] - mesh->xdisp);
	new_vtx->x[1] = mesh->scale * (vertex[1] - mesh->ydisp);
  
	msh_trg_t* trg = mesh_locate_vtx(mesh, new_vtx);

	if (NULL == trg) {
		mvtx_destroy(mesh, new_vtx);
		return false;
	} else {
		nb_ruppert_insert_verified_vtx(mesh, trg, new_vtx);
		return true;
	}
}

void nb_ruppert_insert_verified_vtx(nb_tessellator2D_t *restrict mesh,
				    msh_trg_t *trg,
				    const msh_vtx_t *vtx)
{
	insert_vertex(mesh, trg, vtx, NULL, NULL, NULL);

}

void nb_ruppert_insert_verified_subsgm_midpoint(nb_tessellator2D_t *restrict mesh,
						msh_edge_t *sgm)
{
	msh_vtx_t *v = get_midpoint(mesh, sgm);
	msh_edge_t* subsgm[2];
	insert_midpoint(mesh, sgm, v, NULL, NULL, subsgm, NULL);

}

static void delete_bad_trg(nb_tessellator2D_t *mesh,
			   nb_container_t *encroached_sgm,
			   nb_container_t *big_trg,
			   hash_trg_t *poor_quality_trg)
{
	uint32_t iter = 0;
	while ((nb_container_is_not_empty(big_trg) ||
		hash_trg_length(poor_quality_trg) > 0) &&
	       check_max_trg(mesh) && check_max_vtx(mesh)) {
		/* Re-allocate hash tables if they are too small */
		if (iter % 5000 == 0) {
			if (nb_bins2D_get_min_points_x_bin(mesh->ug_vtx) > 100)
				reallocate_bins(mesh);
		}
		iter ++;

		/* Get next poor-quality triangle */
		msh_trg_t *trg;
		if (nb_container_is_not_empty(big_trg))
			trg = nb_container_delete_first(big_trg);
		else
			trg = hash_trg_remove_first(poor_quality_trg);

		delete_poor_quality_trg(mesh, encroached_sgm, big_trg,
					poor_quality_trg, trg);
		
	}
}

static void delete_poor_quality_trg(nb_tessellator2D_t *mesh,
				    nb_container_t *encroached_sgm,
				    nb_container_t *big_trg,
				    hash_trg_t *poor_quality_trg,
				    msh_trg_t *trg)
{
	/* Get circumcenter */
	msh_vtx_t *cc = mvtx_create(mesh);
	nb_utils2D_get_circumcenter(trg->v1->x, trg->v2->x,
				    trg->v3->x, cc->x);
    
	msh_trg_t *trg_containing_cc =
		get_trg_containing_circumcenter(mesh, trg, cc);

	nb_container_t *sgm_encroached_by_cc =
		nb_allocate_on_stack(nb_container_get_memsize(NB_QUEUE));
	nb_container_init(sgm_encroached_by_cc, NB_QUEUE); 

	get_sgm_encroached_by_vertex(mesh, trg_containing_cc, cc,
				     sgm_encroached_by_cc);

	if (nb_container_is_empty(sgm_encroached_by_cc))
		circumcenter_accepted(mesh, big_trg, poor_quality_trg,
				      trg_containing_cc, cc);
	else
		circumcenter_rejected(mesh, encroached_sgm, big_trg,
				      poor_quality_trg,
				      sgm_encroached_by_cc, cc, trg);

	nb_container_finish(sgm_encroached_by_cc);
}

static void circumcenter_accepted(nb_tessellator2D_t *mesh,
				  nb_container_t *big_trg,
				  hash_trg_t *poor_quality_trg,
				  msh_trg_t *trg_containing_cc,
				  msh_vtx_t *cc)
{
	nb_container_t *new_trg =
		nb_allocate_on_stack(nb_container_get_memsize(NB_QUEUE));
	nb_container_init(new_trg, NB_QUEUE);

	insert_vertex(mesh, trg_containing_cc, cc, big_trg,
		      poor_quality_trg, new_trg);

	verify_new_encroachments(mesh, cc, new_trg, NULL,
				 big_trg, poor_quality_trg);
	nb_container_finish(new_trg);
}

static void circumcenter_rejected(nb_tessellator2D_t *mesh,
				  nb_container_t *encroached_sgm,
				  nb_container_t *big_trg,
				  hash_trg_t *poor_quality_trg,
				  nb_container_t *sgm_encroached_by_cc,
				  msh_vtx_t *cc, msh_trg_t *trg)
{
	mvtx_destroy(mesh, cc);

	/* Process segments encroached by the circumcenter */
	double d = nb_utils2D_get_min_trg_edge(trg->v1->x,
					       trg->v2->x,
					       trg->v3->x);

	while (nb_container_is_not_empty(sgm_encroached_by_cc)) {
		msh_edge_t* restrict sgm =
			nb_container_delete_first(sgm_encroached_by_cc);
		/* Insert segment in encroached list if the triangle is to big 
		 * or if the split is permitted */
		bool is_encroached;
		if (mtrg_is_too_big(mesh, trg, NULL))		
			is_encroached = true;
		else if (split_is_permitted(sgm, d))
			is_encroached = true;
		else
			is_encroached = false;

		if (is_encroached)
			nb_container_insert(encroached_sgm, sgm);
	}
	/* Process encroached segments */
	if (nb_container_is_not_empty(encroached_sgm)) {
		/* Put it back for another try */
		check_trg(trg, mesh, big_trg, poor_quality_trg);
	
		/* Split encroached segments */
		split_encroached_segments(mesh, encroached_sgm,
					  big_trg,
					  poor_quality_trg);
	}
}

static void reallocate_bins(nb_tessellator2D_t *const restrict mesh)
{
	double avg_points_x_bin = 
		nb_bins2D_get_length(mesh->ug_vtx) / 
		nb_bins2D_get_N_bins(mesh->ug_vtx);
	double new_bin_size =
		nb_bins2D_get_size_of_bins(mesh->ug_vtx) /
		sqrt(avg_points_x_bin);
  
	nb_bins2D_t *vertices = nb_allocate_on_stack(nb_bins2D_get_memsize());
	nb_bins2D_init(vertices, new_bin_size);

	while (nb_bins2D_is_not_empty(mesh->ug_vtx)) {
		msh_vtx_t* vtx = nb_bins2D_delete_first(mesh->ug_vtx);
		nb_bins2D_insert(vertices, vtx);
	}
	nb_bins2D_finish(mesh->ug_vtx);
	mesh->ug_vtx = vertices;
}

static inline hash_trg_t* hash_trg_create(void)
{
	hash_trg_t* htrg = nb_allocate_zero_mem(sizeof(hash_trg_t));
	for (uint32_t i = 0; i < 64; i++) {
		htrg->avl[i] = nb_container_create(NB_SORTED);
		nb_container_set_comparer(htrg->avl[i], compare_trg_attr);
	}
	return htrg;
}

static inline uint32_t hash_trg_length(const hash_trg_t *const restrict htrg)
{
	return htrg->length;
}

static inline void hash_trg_insert(hash_trg_t *const restrict htrg,
				   msh_trg_t *const restrict trg)
{
	uint32_t hash_key = ((uint32_t) trg->feature) % 64;

	bool is_inserted =
		nb_container_insert(htrg->avl[hash_key], trg);

	if (is_inserted)
		htrg->length += 1;
}

static inline msh_trg_t* hash_trg_remove_first
                         (hash_trg_t *const restrict htrg)
{
	msh_trg_t *trg = NULL;
	for (int i = 63; i >= 0; i--) {
		if (nb_container_is_not_empty(htrg->avl[i])) {
			htrg->length -= 1;
			trg = nb_container_delete_first(htrg->avl[i]);
			break;
		}
	}
	return trg;
}

static inline void hash_trg_remove(hash_trg_t *const restrict htrg,
				   msh_trg_t *const restrict trg)
{
	uint32_t hash_key = ((uint32_t) trg->feature) % 64;

	bool is_removed = nb_container_delete(htrg->avl[hash_key], trg);

	if (is_removed)
		htrg->length -= 1;
}

static inline void hash_trg_destroy(hash_trg_t* restrict htrg)
{
	for (uint32_t i = 0; i < 64; i++)
		nb_container_destroy(htrg->avl[i]);
	nb_free_mem(htrg);
}

static inline bool check_max_vtx(const nb_tessellator2D_t *const restrict mesh)
{
	bool allow = true;
	if (0 < mesh->max_vtx)
		allow = (nb_bins2D_get_length(mesh->ug_vtx) < mesh->max_vtx);
	return allow;
}

static inline bool check_max_trg(const nb_tessellator2D_t *const restrict mesh)
{
	
	bool allow = true;
	if (0 < mesh->max_trg)
		allow = (nb_container_get_length(mesh->ht_trg) < 
			 mesh->max_trg);
	return allow;
}

static inline bool is_encroached
                  (const msh_edge_t *const restrict sgm)
{
	bool is_encroached = false;
	msh_vtx_t* restrict vtx;
	if (NULL != sgm->t1) {
		vtx = mtrg_get_opposite_vertex_guided(sgm->t1, sgm, true);
		is_encroached = 
			nb_utils2D_pnt_lies_strictly_in_diametral_circle(sgm->v1->x,
									sgm->v2->x,
									vtx->x);
	}

	if (!is_encroached && NULL != sgm->t2) {
		vtx = mtrg_get_opposite_vertex_guided(sgm->t2, sgm, false);
		is_encroached = 
			nb_utils2D_pnt_lies_strictly_in_diametral_circle(sgm->v1->x,
									sgm->v2->x,
									vtx->x);
	}
	return is_encroached;
}

static inline msh_vtx_t* get_midpoint(nb_tessellator2D_t *mesh,
				      const msh_edge_t *const restrict sgm)
{
	/* Calculate the new vertex (using concentric shells) */
	msh_vtx_t *v = mvtx_create(mesh);
	mvtx_set_type_location(v, ONSEGMENT);
  
	/* Use midpoint */
	v->x[0] = 0.5 * (sgm->v1->x[0] + sgm->v2->x[0]);
	v->x[1] = 0.5 * (sgm->v1->x[1] + sgm->v2->x[1]);

	if (mvtx_is_type_origin(sgm->v1, INPUT) ||
	    mvtx_is_type_origin(sgm->v2, INPUT)) 
		concentric_shell(sgm, v);
	return v;
}

static void concentric_shell(const msh_edge_t *const sgm, msh_vtx_t *v)
{
	/* Use the input vertex, stored in v1,  as the centroid of the 
	 * concentric shells */
	msh_vtx_t* v1 = sgm->v1;
	msh_vtx_t* v2 = sgm->v2;
	if (mvtx_is_type_origin(sgm->v1, STEINER)) {
		v1 = sgm->v2;
		v2 = sgm->v1;
	}
	/* Use concentric shells */
	double length = nb_utils2D_get_dist(v->x, v1->x);
	double k_real = nb_math_log2(length / _NB_CC_SHELL_UNIT);
	int k_int = (int)((k_real > 0)?(k_real + 0.5):(k_real - 0.5));
	/* factor = UNIT * 2^k / (2.0 * length) */
	double factor = _NB_CC_SHELL_UNIT * pow(2.0, k_int - 1.0) / length;
	v->x[0] = v1->x[0] * (1.0 - factor) + v2->x[0] * factor;
	v->x[1] = v1->x[1] * (1.0 - factor) + v2->x[1] * factor;
}

static void get_encroached_triangles
                        (const msh_trg_t *const restrict first_trg_to_check,
			 const msh_vtx_t *const restrict v,
			 nb_container_t *encroached_trg)
{
	/* Search encroached triangles */
	nb_container_t *const unencroached_trg =
		nb_allocate_on_stack(nb_container_get_memsize(NB_SORTED));
	nb_container_init(unencroached_trg, NB_SORTED);
	nb_container_t *const processing_trg =
		nb_allocate_on_stack(nb_container_get_memsize(NB_SORTED));
	nb_container_init(processing_trg, NB_SORTED);
  
	nb_container_insert(processing_trg, first_trg_to_check);

	while (nb_container_is_not_empty(processing_trg)) {
		msh_trg_t *trg = nb_container_delete_first(processing_trg);
		/* Check if it is encroached */
		if (nb_utils2D_pnt_lies_strictly_in_circumcircle(trg->v1->x,
								  trg->v2->x,
								  trg->v3->x,
								  v->x)) {
			nb_container_insert(encroached_trg, trg);
			check_encroached_neighbours(trg, encroached_trg,
						    unencroached_trg,
						    processing_trg);
		} else {
			nb_container_insert(unencroached_trg, trg);
		}
	}

	nb_container_finish(unencroached_trg);
	nb_container_finish(processing_trg);
}

static void check_encroached_neighbours(const msh_trg_t *trg,
					nb_container_t *encroached_trg,
					nb_container_t *unencroached_trg,
					nb_container_t *processing_trg)
{
	if (trg->t1 != NULL) {
		if (!medge_is_subsgm(trg->s1))
			if (NOT_EXIST(encroached_trg, trg->t1))
				if (NOT_EXIST(unencroached_trg,	trg->t1))
					nb_container_insert(processing_trg,
							    trg->t1);
	}
	if (trg->t2 != NULL) {
		if (!medge_is_subsgm(trg->s2))
			if (NOT_EXIST(encroached_trg, trg->t2))
				if (NOT_EXIST(unencroached_trg,	trg->t2))
					nb_container_insert(processing_trg,
							    trg->t2);
	}
	if (trg->t3 != NULL){
		if (!medge_is_subsgm(trg->s3))
			if (NOT_EXIST(encroached_trg, trg->t3))
				if (NOT_EXIST(unencroached_trg,	trg->t3))
					nb_container_insert(processing_trg,
							    trg->t3);
	}
}

static void remove_encroached_triangles
             (nb_tessellator2D_t *const restrict mesh,
	      msh_trg_t *const restrict first_trg_to_check,
	      const msh_vtx_t *const restrict v,
	      /* v_orfans_reference: NULL to take an arbitrary reference */
	      const msh_vtx_t *const restrict v_orfans_reference,
	      /* big_trg: NULL if not required */
	      nb_container_t *const restrict big_trg,
	      /* poor_quality_trg: NULL if not required */
	      hash_trg_t *const restrict poor_quality_trg,
	      nb_container_t* orfan_vtx)
{
	nb_container_t *encroached_trg =
		nb_allocate_on_stack(nb_container_get_memsize(NB_SORTED));
	nb_container_init(encroached_trg, NB_SORTED);

	get_encroached_triangles(first_trg_to_check, v, encroached_trg);

	/* Destroy encroached triangles */
	while (nb_container_is_not_empty(encroached_trg)) {
		msh_trg_t *trg = nb_container_delete_first(encroached_trg);
		/* Set vertices as orfans :( */
		nb_container_insert(orfan_vtx, trg->v1);
		nb_container_insert(orfan_vtx, trg->v2);
		nb_container_insert(orfan_vtx, trg->v3);

		mesh_substract_triangle(mesh, trg);

		/* Substract from poor-quality set (if exist) */
		if (NULL != poor_quality_trg)
			hash_trg_remove(poor_quality_trg, trg);

		/* Substract from big triangles set (if exist) */
		if (NULL != big_trg)
			nb_container_delete(big_trg, trg);

		mtrg_nb_free_mem(mesh, trg);
	}
	nb_container_finish(encroached_trg);
}

static inline msh_vtx_t* retrg_fan_get_next_trg
                          (uint32_t* N_orfan_vtx,
			   msh_vtx_t ** orfan_vtx,
			   const msh_vtx_t *const v1,
			   const msh_vtx_t *const v2){
	double reference_angle = atan2(v2->x[1] - v1->x[1],
				       v2->x[0] - v1->x[0]);
	double min_angle = NB_MATH_PI;
	uint32_t i3 = *N_orfan_vtx;
	for (uint32_t i = 0; i < *N_orfan_vtx; i++) {
		if (v2 == orfan_vtx[i])
			continue;
		double angle = atan2(orfan_vtx[i]->x[1] - v1->x[1],
				     orfan_vtx[i]->x[0] - v1->x[0]);
		double angle_diff = angle - reference_angle;
		if (angle_diff < 0.0)
			angle_diff += 2.0 * NB_MATH_PI;
		if (angle_diff < min_angle) {
			i3 = i;
			min_angle = angle_diff;
		}
	}
  
	if (i3 < *N_orfan_vtx) {
		*N_orfan_vtx -= 1;
		msh_vtx_t* aux = orfan_vtx[i3];
		orfan_vtx[i3] = orfan_vtx[*N_orfan_vtx];
		orfan_vtx[*N_orfan_vtx] = aux;
		return aux;
	} else {
		return NULL;
	}
}

static inline void retriangulate_fan(nb_tessellator2D_t *const mesh,
				     const msh_vtx_t *const v_pivot,
				     const msh_vtx_t *const v_start,
				     nb_container_t *orfan_vtx,
				     /* l_new_trg: NULL if not required */
				     nb_container_t *const l_new_trg)
{
	/* Allocate orfan vertices into an array */
	uint32_t orfan_memsize = nb_container_get_length(orfan_vtx)
		* sizeof(msh_vtx_t*);
	msh_vtx_t** orfan_array = nb_soft_allocate_mem(orfan_memsize);

	uint32_t N_orfan_vtx = 0;
	while (nb_container_is_not_empty(orfan_vtx)) {
		msh_vtx_t *vtx = nb_container_delete_first(orfan_vtx);
		orfan_array[N_orfan_vtx] = vtx;
		N_orfan_vtx += 1;
	}

	/* Retriangulate fan */
	msh_vtx_t *vtx = (msh_vtx_t*) v_start;
	msh_edge_t *sgm;
	do {
		msh_vtx_t *restrict v3 = 
			retrg_fan_get_next_trg(&N_orfan_vtx, 
					       orfan_array,
					       v_pivot, vtx);

		/* Create triangle */
		msh_trg_t *trg = mtrg_allocate_zero_mem(mesh);
		trg->v1 = (msh_vtx_t*) v_pivot;
		trg->v2 = vtx;
		trg->v3 = v3;

		/* Add triangle to the mesh */
		mesh_add_triangle(mesh, trg);
		mesh->do_after_insert_trg(mesh);

		if (l_new_trg != NULL)
			nb_container_insert(l_new_trg, trg);

		vtx = v3;

		sgm = mesh_exist_edge(mesh->ht_edge, v_pivot, vtx);
	} while (vtx != v_start && !medge_is_subsgm(sgm));

	nb_soft_free_mem(orfan_memsize, orfan_array);
}

static void verify_new_encroachments
                             (nb_tessellator2D_t *const restrict mesh,
			      const msh_vtx_t *const restrict v,
			      nb_container_t *const restrict l_new_trg,
			      /* Could be NULL if not required */
			      nb_container_t *const restrict encroached_sgm,
			      nb_container_t *const restrict big_trg,
			      hash_trg_t *const restrict poor_quality_trg)
{
	bool (*inside)(const double[2], const double[2], const double [2]) =
		nb_utils2D_pnt_lies_strictly_in_diametral_circle;

	while (nb_container_is_not_empty(l_new_trg)) {
		msh_trg_t* restrict trg = nb_container_delete_first(l_new_trg);
		if (NULL != encroached_sgm) {
			msh_edge_t* restrict edge =
				mtrg_get_opposite_edge(trg, v);
			if (medge_is_subsgm(edge)) {
				if (inside(edge->v1->x, edge->v2->x, v->x)) {
					nb_container_insert(encroached_sgm,
							     edge);
					continue;
				}
			}
		}		
		check_trg(trg, mesh, big_trg, poor_quality_trg);
	}
}

void insert_vertex(nb_tessellator2D_t *mesh,
		   msh_trg_t *const restrict trg_containing_cc,
		   const msh_vtx_t *const restrict cc,
		   /* big_trg: NULL if not required */
		   nb_container_t *const restrict  big_trg,
		   /* poor_quality_trg: NULL if not required */
		   hash_trg_t *const restrict poor_quality_trg,
		   /* l_new_trg: NULL if not required */
		   nb_container_t *const restrict l_new_trg)
{
	nb_bins2D_insert(mesh->ug_vtx, cc);

	nb_container_t* orfan_vtx =
		nb_allocate_on_stack(nb_container_get_memsize(NB_SORTED));
	nb_container_init(orfan_vtx, NB_SORTED);

	remove_encroached_triangles(mesh, trg_containing_cc, 
				    cc, NULL, big_trg,
				    poor_quality_trg,
				    orfan_vtx);
	retriangulate_fan(mesh, cc, 
			  nb_container_get_first(orfan_vtx),
			  orfan_vtx,
			  l_new_trg);
	nb_container_finish(orfan_vtx);
	mesh->do_after_insert_vtx(mesh);
}

static void insert_midpoint(nb_tessellator2D_t *const mesh,
			    msh_edge_t *const sgm,
			    const msh_vtx_t *const v,
			    /* big_trg: NULL if not required */
			    nb_container_t *const big_trg,
			    /* poor_quality_trg: NULL if not required */
			    hash_trg_t *const poor_quality_trg,
			    msh_edge_t* subsgm[2],
			    /* l_new_trg: NULL if not required */
			    nb_container_t *const l_new_trg)
{
	nb_bins2D_insert(mesh->ug_vtx, v);

	/* Split the subsegment */  
	subsgm[0] = medge_allocate_zero_mem(mesh);
	subsgm[0]->v1 = sgm->v1;
	subsgm[0]->v2 = (msh_vtx_t*) v;
	nb_container_insert(mesh->ht_edge, subsgm[0]);
	
	subsgm[1] = medge_allocate_zero_mem(mesh);
	subsgm[1]->v1 = (msh_vtx_t*) v;
	subsgm[1]->v2 = sgm->v2;
	nb_container_insert(mesh->ht_edge, subsgm[1]);

	/* Replace segment by its subsegments */
	uint32_t input_idx = medge_subsgm_idx(sgm);
	msh_edge_t* prev = medge_subsgm_prev(sgm);
	msh_edge_t* next = medge_subsgm_next(sgm);
	medge_set_as_subsgm(subsgm[0], input_idx,
			    prev, subsgm[1]);
	medge_set_as_subsgm(subsgm[1], input_idx,
			    subsgm[0], next);
	if (NULL != prev)
		medge_update_subsgm_next(prev, subsgm[0]);
	if (NULL != next)
		medge_update_subsgm_prev(next, subsgm[1]);
	if (mesh->input_sgm[input_idx] == sgm)
		mesh->input_sgm[input_idx] = subsgm[0];
  
	/* Retriangulate left */
	if (NULL != sgm->t1) {
		nb_container_t* orfan_vtx =
			nb_allocate_on_stack(nb_container_get_memsize(NB_SORTED));
		nb_container_init(orfan_vtx, NB_SORTED);

		remove_encroached_triangles(mesh, sgm->t1, v, sgm->v2,
					    big_trg, poor_quality_trg,
					    orfan_vtx);

		retriangulate_fan(mesh, v, sgm->v2, orfan_vtx, l_new_trg);
		nb_container_finish(orfan_vtx);
	}
	/* Retriangulate right */
	if (NULL != sgm->t2) {
		nb_container_t* orfan_vtx =
			nb_allocate_on_stack(nb_container_get_memsize(NB_SORTED));
		nb_container_init(orfan_vtx, NB_SORTED);

		remove_encroached_triangles(mesh, sgm->t2,
					    v, sgm->v1,
					    big_trg,
					    poor_quality_trg,
					    orfan_vtx);

		retriangulate_fan(mesh, v, sgm->v1, orfan_vtx, l_new_trg);
		nb_container_finish(orfan_vtx);
	}
	/* Free splitted segment */
	nb_container_delete(mesh->ht_edge, sgm);
	medge_destroy_subsgm_attribute(sgm);
	medge_nb_free_mem(mesh, sgm);
}

static void split_encroached_segments
                             (nb_tessellator2D_t *const mesh,
			      nb_container_t *const restrict encroached_sgm,
			      nb_container_t *const restrict big_trg,
			      hash_trg_t *const restrict poor_quality_trg)
{
	nb_container_t* l_new_trg = 
		nb_allocate_on_stack(nb_container_get_memsize(NB_QUEUE));
	while (nb_container_is_not_empty(encroached_sgm)) {
		msh_edge_t* sgm = nb_container_delete_first(encroached_sgm);
		msh_vtx_t *v = get_midpoint(mesh, sgm);
		/* Split the segment */
		nb_container_init(l_new_trg, NB_QUEUE);

		msh_edge_t* subsgm[2];
		insert_midpoint(mesh, sgm, v, big_trg, poor_quality_trg,
				subsgm, l_new_trg);

		/* Verify new encroachments */
		verify_new_encroachments(mesh, v, l_new_trg, encroached_sgm,
					 big_trg, poor_quality_trg);
		nb_container_finish(l_new_trg);

		/* Check if the first subsegment is encroached */
		if (is_encroached(subsgm[0]))
			nb_container_insert(encroached_sgm, subsgm[0]);

		/* Check if the second subsegment is encroached */
		if (is_encroached(subsgm[1]))
			nb_container_insert(encroached_sgm, subsgm[1]);

		mesh->do_after_insert_vtx(mesh);
	}
}

static inline void initialize_encroached_sgm
                          (nb_tessellator2D_t *const restrict mesh,
			   nb_container_t *const restrict encroached_sgm)
{
	for (uint32_t i = 0; i < mesh->N_input_sgm; i++) {
		msh_edge_t* restrict sgm = mesh->input_sgm[i];
		while (sgm != NULL) {
			if (is_encroached(sgm))
				nb_container_insert(encroached_sgm, sgm);
			sgm = medge_subsgm_next(sgm);
		}
	}
}

static inline void check_trg(msh_trg_t *const restrict trg,
			     nb_tessellator2D_t *const restrict mesh,
			     nb_container_t *const restrict big_trg,
			     hash_trg_t *const restrict poor_quality_trg)
{
	double big_ratio;
	if (mtrg_is_too_big(mesh, trg, &big_ratio)) {
		trg->feature = big_ratio;
		nb_container_insert(big_trg, trg);
	} else {
		if (1e30 > mesh->cr2se_ratio) {
			double cr2se_ratio =
				nb_utils2D_get_cr2se_ratio(trg->v1->x,
							    trg->v2->x,
							    trg->v3->x);
			if (cr2se_ratio > mesh->cr2se_ratio) {
				trg->feature = cr2se_ratio;
				hash_trg_insert(poor_quality_trg, trg);
			}
		}
	}
}

static void initialize_big_and_poor_quality_trg
                          (nb_tessellator2D_t *const restrict mesh,
			   nb_container_t *const restrict big_trg,
			   hash_trg_t *const restrict poor_quality_trg)
{
	uint32_t size = nb_iterator_get_memsize();
	nb_iterator_t *iter = nb_allocate_on_stack(size);
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, mesh->ht_trg);
	while (nb_iterator_has_more(iter)) {
		msh_trg_t *restrict trg =
			(msh_trg_t*) nb_iterator_get_next(iter);
		check_trg(trg, mesh, big_trg, poor_quality_trg);
	}
	nb_iterator_finish(iter);
}

static inline msh_trg_t* get_trg_containing_circumcenter
                          (const nb_tessellator2D_t *const restrict mesh,
			   const msh_trg_t *const restrict trg,
			   const msh_vtx_t *const restrict cc)
{
	/* Walk direct to the circumcenter */
	msh_trg_t *restrict prev2_trg = NULL;
	msh_trg_t *restrict prev1_trg = NULL;
	msh_trg_t *restrict trg_containing_cc = (msh_trg_t*) trg;
	while (!nb_utils2D_pnt_lies_in_trg(trg_containing_cc->v1->x,
					    trg_containing_cc->v2->x,
					    trg_containing_cc->v3->x,
					    cc->x)) {
		/* Update history */
		prev2_trg = prev1_trg;
		prev1_trg = trg_containing_cc;

		/* Calculate centroid of the triangle */
		double centroid[2];
		nb_utils2D_trg_get_centroid(trg_containing_cc->v1->x,
					     trg_containing_cc->v2->x,
					     trg_containing_cc->v3->x,
					     centroid);

		/* Move to the next triangle */
		if (NB_INTERSECTED == nb_utils2D_get_sgm_intersection
		    				(centroid, cc->x,
						 trg_containing_cc->v1->x,
						 trg_containing_cc->v2->x,
						 NULL))
			trg_containing_cc = trg_containing_cc->t1;
		else if (NB_INTERSECTED == nb_utils2D_get_sgm_intersection
						(centroid, cc->x,
						 trg_containing_cc->v2->x,
						 trg_containing_cc->v3->x,
						 NULL))
			trg_containing_cc = trg_containing_cc->t2;
		else
			trg_containing_cc = trg_containing_cc->t3;

		/* Handle numerical errors: When the cc lies on the boundary
		 * between triangles or outside the mesh. 
		 * The first condition is to prevent a eternal iteration. */
		if (trg_containing_cc == prev2_trg || trg_containing_cc == NULL)
			break;
	}
	return trg_containing_cc;
}

static void get_sgm_encroached_by_vertex
                          (const nb_tessellator2D_t *const restrict mesh,
			   const msh_trg_t *const restrict trg_containing_vtx,
			   const msh_vtx_t *const restrict vtx,
			   nb_container_t *encroached_sgm)
{
	/* Get encroached triangles */
	nb_container_t *encroached_trg =
		nb_allocate_on_stack(nb_container_get_memsize(NB_SORTED));
	nb_container_init(encroached_trg, NB_SORTED);

	get_encroached_triangles(trg_containing_vtx, vtx, encroached_trg);

	/* Get segments */
	nb_container_t *segments =
		nb_allocate_on_stack(nb_container_get_memsize(NB_SORTED));
	nb_container_init(segments, NB_SORTED);

	while (nb_container_is_not_empty(encroached_trg)) {
		msh_trg_t *const restrict trg = 
			nb_container_delete_first(encroached_trg);
		if (medge_is_subsgm(trg->s1))
			nb_container_insert(segments, trg->s1);
		if (medge_is_subsgm(trg->s2))
			nb_container_insert(segments, trg->s2);
		if (medge_is_subsgm(trg->s3))
			nb_container_insert(segments, trg->s3);
	}
	/* Free memory */
	nb_container_finish(encroached_trg);

	/* Get encroached segments */
	while (nb_container_is_not_empty(segments)) {
		msh_edge_t *const restrict sgm = 
			nb_container_delete_first(segments);
		if (nb_utils2D_pnt_lies_strictly_in_diametral_circle(sgm->v1->x,
								      sgm->v2->x,
								      vtx->x))
			nb_container_insert(encroached_sgm, sgm);
	}

	/* Free memory */
	nb_container_finish(segments);
}

static inline bool split_is_permitted(const msh_edge_t *const restrict sgm,
				      double d)
{
	/* Check if the length is not a power of two */
	double length = nb_utils2D_get_dist(sgm->v1->x, sgm->v2->x);
	double k_real = nb_math_log2(length / _NB_CC_SHELL_UNIT);
	int k_int = (int)((k_real > 0)?(k_real + 0.5):(k_real - 0.5));
	double shell_length = _NB_CC_SHELL_UNIT * pow(2.0, k_int);
	if (fabs(length - shell_length) > NB_GEOMETRIC_TOL)
		return true;

	/* Get subsegment clusters */
	double smallest_angle;
	nb_container_t* cluster1 =
		nb_allocate_on_stack(nb_container_get_memsize(NB_SORTED));
	nb_container_init(cluster1, NB_SORTED);

	get_subsgm_cluster(sgm, sgm->v1, &smallest_angle, cluster1);

	double smallest_angle_aux;
	nb_container_t* cluster2 =
		nb_allocate_on_stack(nb_container_get_memsize(NB_SORTED));
	nb_container_init(cluster2, NB_SORTED);

	get_subsgm_cluster(sgm, sgm->v2, &smallest_angle_aux, cluster2);

	/* Check if s belongs to no cluster or if sgm belongs to 
	 * both clusters */
	if (nb_container_is_empty(cluster1) == 
	    nb_container_is_empty(cluster2)) {
		nb_container_finish(cluster1);
		nb_container_finish(cluster2);
		return true;
	}

	/* Get the cluster containing the sgm */
	nb_container_t *cluster;
	if (nb_container_is_not_empty(cluster1)) {
		cluster = cluster1;
		nb_container_finish(cluster2);
	} else {
		cluster = cluster2;
		nb_container_finish(cluster1);
		smallest_angle = smallest_angle_aux;
	}

	/* Check if the cluster contains a segment shorter than sgm */
	while (nb_container_is_not_empty(cluster)) {
		const msh_edge_t *const restrict cluster_sgm = 
			nb_container_delete_first(cluster);

		if (cluster_sgm == sgm) 
			continue;
		double length_aux = nb_utils2D_get_dist(cluster_sgm->v1->x,
					     cluster_sgm->v2->x);
		if (length - length_aux > NB_GEOMETRIC_TOL) {
			nb_container_finish(cluster);
			return true;
		}
	}
	nb_container_finish(cluster);
  
	/* No new edge will be shorter than shortest existing edge */
	double r_min = length * sin(smallest_angle / 2.0);
	return r_min >= d + NB_GEOMETRIC_TOL;
}

static void get_subsgm_cluster(const msh_edge_t *const sgm,
			       const msh_vtx_t *const vtx,
			       double* smallest_angle,
			       nb_container_t* cluster)
{
	*smallest_angle = NB_MATH_PI;

	/* Get subsegments from the left */
	msh_vtx_t* v2 = (sgm->v1 == vtx)? sgm->v2: sgm->v1;
	double angle_prev = atan2(v2->x[1] - vtx->x[1],
			    v2->x[0] - vtx->x[0]);
	msh_edge_t* restrict sgm_next = 
		medge_get_CCW_subsgm(sgm, vtx);
	while (sgm_next != sgm && sgm_next != NULL) {
		v2 = (sgm_next->v1 == vtx)? sgm_next->v2: sgm_next->v1;
		double angle_next = atan2(v2->x[1] - vtx->x[1],
					  v2->x[0] - vtx->x[0]);
		/* Calculate angle between segments */
		double angle = angle_next - angle_prev;
		if (angle < 0.0)
			angle += 2.0 * NB_MATH_PI;
    
		/* If the angle is lower than 60 degrees, add to cluster */
		if (angle < NB_MATH_PI / 3.0) {
			nb_container_insert(cluster, sgm_next);
			/* Get smallest angle */
			if (angle < *smallest_angle)
				*smallest_angle = angle;
    
			/* Get next edge in the fan */
			angle_prev = angle_next;
			sgm_next = medge_get_CCW_subsgm(sgm_next, vtx);
		} else {
			break;
		}
	}
  
	/* Get subsegments from the right */
	v2 = (sgm->v1 == vtx)? sgm->v2: sgm->v1;
	angle_prev = atan2(v2->x[1] - vtx->x[1],
			   v2->x[0] - vtx->x[0]);
	sgm_next = medge_get_CW_subsgm(sgm, vtx);
	while (sgm_next != sgm && sgm_next != NULL) {
		v2 = (sgm_next->v1 == vtx)? sgm_next->v2: sgm_next->v1;
		double angle_next = atan2(v2->x[1] - vtx->x[1],
					  v2->x[0] - vtx->x[0]);
		/* Calculate angle between segments */
		double angle = angle_prev - angle_next;
		if(angle < 0.0)
			angle += 2.0 * NB_MATH_PI;
    
		/* If the angle is lower than 60 degrees, add to cluster */
		if (angle < NB_MATH_PI / 3.0) {
			nb_container_insert(cluster, sgm_next);
			/* Get smallest angle */
			if (angle < *smallest_angle)
				*smallest_angle = angle;
    
			/* Get next edge in the fan */
			angle_prev = angle_next;
			sgm_next = medge_get_CW_subsgm(sgm_next, vtx);
		} else {
			break;
		}
	}
    
	/* Return cluster */
	if (nb_container_is_not_empty(cluster))
		nb_container_insert(cluster, sgm);
	else
		*smallest_angle = 0.0;
}

static bool edge_violates_constrain(const nb_tessellator2D_t *const restrict mesh,
				    const msh_edge_t *const restrict sgm,
				    double *big_ratio)
{
	bool violates_constrain = false;
	double d = nb_utils2D_get_dist(sgm->v1->x, sgm->v2->x);
	*big_ratio = d;
	double edge_max = mesh->max_edge_length * mesh->scale;
	if (d > NB_GEOMETRIC_TOL && edge_max < d &&
	    mesh->has_edge_length_constraint) {
		violates_constrain = true;
	} else if (medge_is_subsgm(sgm) && mesh->has_subsgm_length_constraint) {
		double sgm_max = mesh->max_subsgm_length * mesh->scale;
		if (d > NB_GEOMETRIC_TOL && sgm_max < d) {
			violates_constrain = true;
		}
	}
	return violates_constrain;
}

static bool edge_greater_than_density(const nb_tessellator2D_t *const restrict mesh,
				      const msh_edge_t *const restrict sgm,
				      double *big_ratio)
/* Calculate adimensional length:
 *               
 *            1 _             x(t) = v1 + t * (v2-v1)
 *             |   1          e.g.
 *    lh = |s| | ----- dt,          x( 0 ) = v1
 *            _|  h(x)              x(0.5) = 0.5 * (v1 + v2)
 *            0                     x( 1 ) = v2
 *
 *    Where |s| is the length of the segment, v1 and v2 are
 *    the vertices forming such a segment, and h(·) = 1/density(·)
 *
 *    All the edges of the utopic perfect mesh should have lh = 1.
 */
{
	bool is_greater;
	double d = nb_utils2D_get_dist(sgm->v1->x, sgm->v2->x);
	if (NB_GEOMETRIC_TOL > d) {
		/* Prevent of eternal running */
		is_greater = false;  
	} else {
		d /= mesh->scale;
		double lh = calculate_lh(mesh, sgm, d, 1);
		if (_NB_MAX_LH_TOLERATED < lh) {
			is_greater = true;
		} else {
			lh = calculate_lh(mesh, sgm, d, 2);
			is_greater = (_NB_MAX_LH_TOLERATED < lh);
		}
		*big_ratio = lh;
	}
	return is_greater;
}

static double calculate_lh(const nb_tessellator2D_t *const restrict mesh,
			   const msh_edge_t *const restrict sgm,
			   double dist, int N_trapezoids)
 /*    For convenience, we estimate using the trapezoidal rule
 *
 *                   1 _             x(t) = v1 + t * (v2-v1)
 *                /  |              e.g.
 *    lh =  |s|  /   |  h(x) dt,          x( 0 ) = v1
 *              /   _|                    x(0.5) = 0.5 * (v1 + v2)
 *                  0                     x( 1 ) = v2
 */
{
	double width = 1.0 / N_trapezoids;
	double integral = 0.0;
	double h1 = calculate_h(mesh, sgm, 0);
	for (int i = 0; i < N_trapezoids; i++) {
		double t = (i + 1) * width;
		double h2 = calculate_h(mesh, sgm, t);
		if (h1 > _NB_MAX_GRADING_RATIO * h2)
			h1 = _NB_MAX_GRADING_RATIO * h2;
		else if (h2 > _NB_MAX_GRADING_RATIO * h1)
			h2 = _NB_MAX_GRADING_RATIO * h1;
		
		integral += (h1 + h2);
		h1 = h2;
	}
	return 2 * dist / (integral * width + NB_GEOMETRIC_TOL);
}

static inline double calculate_h(const nb_tessellator2D_t *const restrict mesh,
				 const msh_edge_t *const restrict sgm,
				 double t)
{
	double vtx[2];
	vtx[0] = (1 - t) * sgm->v1->x[0] + t * sgm->v2->x[0];
	vtx[1] = (1 - t) * sgm->v1->x[1] + t * sgm->v2->x[1];
	double vtx_scaled[2];
	mesh_get_extern_scale_and_disp(mesh, vtx, vtx_scaled);
	double density = mesh->density(vtx_scaled, mesh->density_data);
	return 1.0 / (density + NB_GEOMETRIC_TOL);
}

static inline bool mtrg_is_too_big(const nb_tessellator2D_t *const restrict mesh,
				   const msh_trg_t *const restrict trg,
				   /* big_ratio could be NULL if not required */
				   double *big_ratio)
{
	double size_ratio = 1.0;
	bool is_too_big = false;
	if (mesh->has_edge_length_constraint ||
	    mesh->has_subsgm_length_constraint) {
		is_too_big = edge_violates_constrain(mesh, trg->s1, &size_ratio);
		if (!is_too_big)
			is_too_big = edge_violates_constrain(mesh, trg->s2,
							     &size_ratio);
		if (!is_too_big)
			is_too_big = edge_violates_constrain(mesh, trg->s3,
							     &size_ratio);
	}
	if (!is_too_big) {
		if (NULL != mesh->density) {
			is_too_big = edge_greater_than_density(mesh, trg->s1,
							       &size_ratio);
			if (!is_too_big)
				is_too_big = 
					edge_greater_than_density(mesh,
								  trg->s2,
								  &size_ratio);
			if (!is_too_big)
				is_too_big = 
					edge_greater_than_density(mesh,
								  trg->s3,
								  &size_ratio);
		}
	}
	if(NULL != big_ratio)
		*big_ratio = size_ratio;
	return is_too_big;
}
