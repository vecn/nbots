#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "vcn/math_bot.h"
#include "vcn/container_bot.h"
#include "vcn/geometric_bot/utils2D.h"
#include "vcn/geometric_bot/mesh/mesh2D.h"

#include "../mesh2D_structs.h"
#include "ruppert.h"

#define _VCN_MAX_LH_TOLERATED (1.0)
#define _VCN_MAX_GRADING_RATIO (27.0)
#define _VCN_SUBSEGMENT_VTX ((void*)0x2)
#define _VCN_CC_SHELL_UNIT (1e-3)
typedef struct {
	vcn_container_t* avl[64];
	uint32_t length;
} hash_trg_t;

static inline uint32_t key_trg_attr(const void *const trg_ptr);
static hash_trg_t* hash_trg_create(void);
static uint32_t hash_trg_length(const hash_trg_t *const htrg);
static void hash_trg_insert(hash_trg_t *const htrg,
				   msh_trg_t *const trg,
				   double cr2se_ratio);
static msh_trg_t* hash_trg_remove_first(hash_trg_t *const htrg);
static void hash_trg_remove(hash_trg_t *const htrg,
			    msh_trg_t *const trg);
static void hash_trg_destroy(hash_trg_t *const htrg);

static void reallocate_bins(vcn_mesh_t *const restrict mesh);
static bool check_max_vtx(const vcn_mesh_t *const mesh);
static bool check_max_trg(const vcn_mesh_t *const mesh);
static bool is_encroached(const msh_edge_t *const sgm);
static msh_vtx_t* get_midpoint(const msh_edge_t *const sgm);
static vcn_container_t* get_encroached_triangles
                             (const msh_trg_t *const first_trg_to_check,
			      const msh_vtx_t *const v);
static vcn_container_t* remove_encroached_triangles
             (vcn_mesh_t *const mesh,
	      msh_trg_t *const first_trg_to_check,
	      const msh_vtx_t *const v,
	      /* v_orfans_reference: NULL to take an arbitrary reference */
	      const msh_vtx_t *const v_orfans_reference,
	      /* big_trg: NULL if not required */
	      vcn_container_t *const big_trg,
	      /* poor_quality_trg: NULL if not required */
	      hash_trg_t *const poor_quality_trg);

static msh_vtx_t* retrg_fan_get_next_trg
                          (uint32_t* N_orfan_vtx,
			   msh_vtx_t ** orfan_vtx,
			   const msh_vtx_t *const v1,
			   const msh_vtx_t *const v2);

static void retriangulate_fan
                             (vcn_mesh_t *const mesh,
			      const msh_vtx_t *const v_pivot,
			      const msh_vtx_t *const v_start,
			      vcn_container_t *orfan_vtx,
			      /* l_new_trg: NULL if not required */
			      vcn_container_t *const l_new_trg);
static void insert_vertex(vcn_mesh_t *const mesh,
				  msh_trg_t *const trg_containing_cc,
			  const msh_vtx_t *const cc,
			  /* big_trg: NULL if not required */
			  vcn_container_t *const big_trg,
			  /* poor_quality_trg: NULL if not required */
			  hash_trg_t *const poor_quality_trg,
			  /* l_new_trg: NULL if not required */
			  vcn_container_t *const l_new_trg);

static void insert_midpoint(vcn_mesh_t *const mesh,
			    msh_edge_t *const sgm,
			    const msh_vtx_t *const v,
			    vcn_container_t *const big_trg,
			    hash_trg_t *const poor_quality_trg,
			    msh_edge_t* subsgm[2],
			    vcn_container_t *const l_new_trg);

static void split_encroached_segments(vcn_mesh_t *const mesh,
				      vcn_container_t *const encroached_sgm,
				      vcn_container_t *const big_trg,
				      hash_trg_t *const poor_quality_trg);
						   
static void verify_new_encroachments
                             (vcn_mesh_t *const mesh,
			      const msh_vtx_t *const  v,
			      vcn_container_t *const l_new_trg,
			      vcn_container_t *const encroached_sgm,
			      vcn_container_t *const big_trg,
			      hash_trg_t *const poor_quality_trg);

static void initialize_encroached_sgm
                          (vcn_mesh_t *const mesh,
			   vcn_container_t *const encroached_sgm);

static void insert_big_trg(vcn_container_t *big_trg, msh_trg_t *trg,
			   double big_ratio);

static msh_trg_t* remove_bigger_trg(vcn_container_t *big_trg);

static void remove_big_trg(vcn_container_t *big_trg, msh_trg_t *trg);

static void check_trg(msh_trg_t *const trg, vcn_mesh_t *const mesh,
			   vcn_container_t *const big_trg,
			   hash_trg_t *const poor_quality_trg);

static void initialize_big_and_poor_quality_trg
                          (vcn_mesh_t *const mesh,
			   vcn_container_t *const big_trg,
			   hash_trg_t *const poor_quality_trg);

static msh_trg_t* get_trg_containing_circumcenter
                          (const vcn_mesh_t *const mesh,
			   const msh_trg_t *const trg,
			   const msh_vtx_t *const cc);

static vcn_container_t* get_sgm_encroached_by_vertex
                          (const vcn_mesh_t *const mesh,
			   const msh_trg_t *const trg_containing_vtx,
			   const msh_vtx_t *const vtx);

static bool split_is_permitted
                          (const msh_edge_t *const sgm,
			   double d);

static vcn_container_t* get_subsgm_cluster 
                          (const msh_edge_t *const sgm,
			   const msh_vtx_t *const vtx,
			   double* smallest_angle);

static void delete_bad_trg(vcn_mesh_t *mesh,
			   vcn_container_t *encroached_sgm,
			   vcn_container_t *big_trg,
			   hash_trg_t *poor_quality_trg);

static bool has_edge_length_constrained(const vcn_mesh_t *const mesh);
static bool edge_violates_constrain(const vcn_mesh_t *const restrict mesh,
				    const msh_edge_t *const restrict sgm,
				    /* big_ratio can be NULL if not required */
				    double *big_ratio);
static bool edge_greater_than_density(const vcn_mesh_t *const restrict mesh,
				      const msh_edge_t *const restrict sgm,
				      /* big_ratio can be NULL if not required */
				      double *big_ratio);
static double calculate_lh(const vcn_mesh_t *const restrict mesh,
			   const msh_edge_t *const restrict sgm,
			   double dist, int N_trapezoids);
static double calculate_h(const vcn_mesh_t *const restrict mesh,
			  const msh_edge_t *const restrict sgm,
			  double t);
static bool mtrg_is_too_big(const vcn_mesh_t *const restrict mesh,
			    const msh_trg_t *const restrict trg,
			    /* big_ratio could be NULL if not required */
			    double *big_ratio);

void vcn_ruppert_refine(vcn_mesh_t *restrict mesh)
{
	/* Allocate data structures to allocate encroached elements */
	vcn_container_t *restrict encroached_sgm =
		vcn_container_create(VCN_CONTAINER_SORTED);
	vcn_container_set_key_generator(encroached_sgm, hash_key_edge);

	vcn_container_t *restrict big_trg =
		vcn_container_create(VCN_CONTAINER_SORTED);
	vcn_container_set_key_generator(big_trg, key_trg_attr);

	hash_trg_t *restrict poor_quality_trg = hash_trg_create();

	/* Initialize FIFO with encroached segments */
	initialize_encroached_sgm(mesh, encroached_sgm);

	/* Calculate max circumradius to shortest edge ratio allowed */
	if (mesh->min_angle > VCN_MESH_MAX_ANGLE) {
		if (0 == mesh->max_vtx) {
			printf("WARNING in vcn_mesh_refine(): ");
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
	vcn_container_destroy(encroached_sgm);
	vcn_container_destroy(big_trg);
	hash_trg_destroy(poor_quality_trg);
}

static inline uint32_t key_trg_attr(const void *const trg_ptr)
{
	const msh_trg_t *const trg = trg_ptr;
	return *(uint32_t*)(trg->attr);
}

bool vcn_ruppert_insert_vtx(vcn_mesh_t *restrict mesh, const double vertex[2])
{
	msh_vtx_t* new_vtx = (msh_vtx_t*) vcn_point2D_create();

	/* Move and scale new vertex */
	new_vtx->x[0] = mesh->scale * (vertex[0] - mesh->xdisp);
	new_vtx->x[1] = mesh->scale * (vertex[1] - mesh->ydisp);
  
	msh_trg_t* trg = mesh_locate_vtx(mesh, new_vtx);

	if (NULL == trg) {
		vcn_point2D_destroy(new_vtx);
		return false;
	}

	insert_vertex(mesh, trg, new_vtx, NULL, NULL, NULL);
	return true;
}

static void delete_bad_trg(vcn_mesh_t *mesh,
			   vcn_container_t *encroached_sgm,
			   vcn_container_t *big_trg,
			   hash_trg_t *poor_quality_trg)
{
	/* Remove poor quality triangles */
	uint32_t iter = 0;
	while ((vcn_container_is_not_empty(big_trg) ||
		hash_trg_length(poor_quality_trg) > 0) &&
	       check_max_trg(mesh) && check_max_vtx(mesh)) {
		/* Re-allocate hash tables if they are too small */
		if (iter % 5000 == 0) {
			if (vcn_bins2D_get_min_points_x_bin(mesh->ug_vtx) > 100)
				reallocate_bins(mesh);
		}
		iter ++;
		/* Get next poor-quality triangle */
		msh_trg_t *trg;
		if (vcn_container_is_not_empty(big_trg))
			trg = remove_bigger_trg(big_trg);
		else
			trg = hash_trg_remove_first(poor_quality_trg);

		/* Get circumcenter */
		msh_vtx_t *restrict cc = (msh_vtx_t*) vcn_point2D_create();
		vcn_utils2D_get_circumcenter(trg->v1->x, trg->v2->x, trg->v3->x, cc->x);
    
		msh_trg_t *restrict trg_containing_cc =
			get_trg_containing_circumcenter(mesh, trg, cc);

		vcn_container_t *restrict  sgm_encroached_by_cc = 
			get_sgm_encroached_by_vertex(mesh, trg_containing_cc, cc);
		if (vcn_container_is_empty(sgm_encroached_by_cc)) {
			/* Circumcenter accepted */
			vcn_container_t *restrict l_new_trg = 
				vcn_container_create(VCN_CONTAINER_QUEUE);

			insert_vertex(mesh, trg_containing_cc, cc, big_trg,
				      poor_quality_trg, l_new_trg);

			verify_new_encroachments(mesh, cc, l_new_trg, NULL,
						 big_trg, poor_quality_trg);
			vcn_container_destroy(l_new_trg);
		} else {
			/* Circumcenter rejected */
			free(cc);

			/* Process segments encroached by the circumcenter */
			double d = vcn_utils2D_get_min_trg_edge(trg->v1->x,
								trg->v2->x,
								trg->v3->x);

			while (vcn_container_is_not_empty(sgm_encroached_by_cc)) {
				msh_edge_t* restrict sgm =
					vcn_container_delete_first(sgm_encroached_by_cc);
				/* Insert segment in encroached list if the triangle is to big 
				 * or if the split is permitted */
				if (mtrg_is_too_big(mesh, trg, NULL))
					vcn_container_insert(encroached_sgm, sgm);
				else if (split_is_permitted(sgm, d))
					vcn_container_insert(encroached_sgm, sgm);
			}
			/* Process encroached segments */
			if (vcn_container_is_not_empty(encroached_sgm)) {
				/* Put it back for another try */
				check_trg(trg, mesh, big_trg, poor_quality_trg);
	
				/* Split encroached segments */
				split_encroached_segments(mesh, encroached_sgm,
							  big_trg,
							  poor_quality_trg);
			}
		}
		vcn_container_destroy(sgm_encroached_by_cc);
	}
}

static void reallocate_bins(vcn_mesh_t *const restrict mesh)
{
	double avg_points_x_bin = 
		vcn_bins2D_get_length(mesh->ug_vtx) / 
		vcn_bins2D_get_N_bins(mesh->ug_vtx);
	double new_bin_size =
		vcn_bins2D_get_size_of_bins(mesh->ug_vtx) /
		sqrt(avg_points_x_bin);
  
	vcn_bins2D_t *vertices = vcn_bins2D_create(new_bin_size);
	while (vcn_bins2D_is_not_empty(mesh->ug_vtx)) {
		msh_vtx_t* vtx = vcn_bins2D_delete_first(mesh->ug_vtx);
		vcn_bins2D_insert(vertices, vtx);
	}
	vcn_bins2D_destroy(mesh->ug_vtx);
	mesh->ug_vtx = vertices;
}
static inline hash_trg_t* hash_trg_create(void)
{
	hash_trg_t* htrg = calloc(1, sizeof(hash_trg_t));
	for (uint32_t i = 0; i < 64; i++) {
		htrg->avl[i] = vcn_container_create(VCN_CONTAINER_SORTED);
		vcn_container_set_key_generator(htrg->avl[i], key_trg_attr);
	}
	return htrg;
}

static inline uint32_t hash_trg_length(const hash_trg_t *const restrict htrg)
{
	return htrg->length;
}

static inline void hash_trg_insert(hash_trg_t *const restrict htrg,
				   msh_trg_t *const restrict trg,
				   double cr2se_ratio)
{
	if (NULL == trg->attr) {
		/* It is not inserted */
		uint32_t hash_key = (uint32_t) cr2se_ratio;
		if (hash_key > 63)
			hash_key = 63;
  
		/* Set circumradius to shortest edge ratio as attribute */
		uint32_t* attr = malloc(sizeof(uint32_t));
		attr[0] = (uint32_t)(cr2se_ratio * 1e6);
		trg->attr = attr;

		bool is_inserted =
			vcn_container_insert(htrg->avl[hash_key], trg);

		if (is_inserted) {
			htrg->length += 1;
		} else {
			free(trg->attr);
			trg->attr = NULL;
		}
	}
}

static inline msh_trg_t* hash_trg_remove_first
                         (hash_trg_t *const restrict htrg)
{
	for (int i = 63; i >= 0; i--) {
		if (vcn_container_is_not_empty(htrg->avl[i])) {
			htrg->length -= 1;
			msh_trg_t* restrict trg = 
				vcn_container_delete_first(htrg->avl[i]);
			free(trg->attr);
			trg->attr = NULL;
			return trg;
		}
	}
	return NULL;
}

static inline void hash_trg_remove(hash_trg_t *const restrict htrg,
				   msh_trg_t *const restrict trg)
{
	if (NULL != trg->attr) {
		double cr2se_ratio = *((uint32_t*)trg->attr) / 1e6;
		uint32_t hash_key = (uint32_t) cr2se_ratio;
		if (hash_key > 63)
			hash_key = 63;

		bool is_removed = vcn_container_delete(htrg->avl[hash_key], trg);

		if (is_removed) {
			htrg->length -= 1;
			free(trg->attr);
			trg->attr = NULL;
		}
	}
}

static inline void hash_trg_destroy(hash_trg_t* restrict htrg)
{
  for (uint32_t i = 0; i < 64; i++) {
    while (vcn_container_is_not_empty(htrg->avl[i])) {
      msh_trg_t* trg = vcn_container_delete_first(htrg->avl[i]);
      free(trg->attr);
      trg->attr = NULL;      
    }
    vcn_container_destroy(htrg->avl[i]);
  }
  free(htrg);
}

static inline bool check_max_vtx(const vcn_mesh_t *const restrict mesh)
{
	bool allow = true;
	if (0 < mesh->max_vtx)
		allow = (vcn_bins2D_get_length(mesh->ug_vtx) < mesh->max_vtx);
	return allow;
}

static inline bool check_max_trg(const vcn_mesh_t *const restrict mesh)
{
	
	bool allow = true;
	if (0 < mesh->max_trg)
		allow = (vcn_container_get_length(mesh->ht_trg) < 
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
			vcn_utils2D_pnt_lies_strictly_in_diametral_circle(sgm->v1->x,
									sgm->v2->x,
									vtx->x);
	}

	if (!is_encroached && NULL != sgm->t2) {
		vtx = mtrg_get_opposite_vertex_guided(sgm->t2, sgm, false);
		is_encroached = 
			vcn_utils2D_pnt_lies_strictly_in_diametral_circle(sgm->v1->x,
									sgm->v2->x,
									vtx->x);
	}
	return is_encroached;
}

static inline msh_vtx_t* get_midpoint
                         (const msh_edge_t *const restrict sgm)
{
	/* Calculate the new vertex (using concentric shells) */
	msh_vtx_t *v = (msh_vtx_t*) vcn_point2D_create();
	v->attr = _VCN_SUBSEGMENT_VTX; 
  
	/* Use midpoint */
	v->x[0] = 0.5 * (sgm->v1->x[0] + sgm->v2->x[0]);
	v->x[1] = 0.5 * (sgm->v1->x[1] + sgm->v2->x[1]);

	if (sgm->v1->attr != _VCN_INPUT_VTX && sgm->v2->attr != _VCN_INPUT_VTX)
		return v;
  
	/* Use the input vertex, stored in v1,  as the centroid of the 
	 * concentric shells */
	msh_vtx_t* v1 = sgm->v1;
	msh_vtx_t* v2 = sgm->v2;
	if (sgm->v1->attr != _VCN_INPUT_VTX) {
		v1 = sgm->v2;
		v2 = sgm->v1;
	}

	/* Use concentric shells */
	double length = vcn_utils2D_get_dist(v->x, v1->x);
	double k_real = vcn_math_log2(length / _VCN_CC_SHELL_UNIT);
	int k_int = (int)((k_real > 0)?(k_real + 0.5):(k_real - 0.5));
	/* factor = UNIT * 2^k / (2.0 * length) */
	double factor = _VCN_CC_SHELL_UNIT * pow(2.0, k_int - 1.0) / length;
	v->x[0] = v1->x[0] * (1.0 - factor) + v2->x[0] * factor;
	v->x[1] = v1->x[1] * (1.0 - factor) + v2->x[1] * factor;
	return v;
}

static inline vcn_container_t* get_encroached_triangles
                        (const msh_trg_t *const restrict first_trg_to_check,
			 const msh_vtx_t *const restrict v)
{
	/* Search encroached triangles */
	vcn_container_t *encroached_trg =
		vcn_container_create(VCN_CONTAINER_SORTED);
	vcn_container_t *const restrict unencroached_trg =
		vcn_container_create(VCN_CONTAINER_SORTED);
	vcn_container_t *const restrict processing_trg =
		vcn_container_create(VCN_CONTAINER_SORTED);
  
	vcn_container_insert(processing_trg, first_trg_to_check);

	while (vcn_container_is_not_empty(processing_trg)) {
		msh_trg_t *const restrict trg = 
			vcn_container_delete_first(processing_trg);
		/* Check if it is encroached */
		if (vcn_utils2D_pnt_lies_strictly_in_circumcircle(trg->v1->x,
								  trg->v2->x,
								  trg->v3->x,
								  v->x)) {
			vcn_container_insert(encroached_trg, trg);
			/* Set the neighbouring triangles to be processed */
			if (trg->t1 != NULL) {
				if (!medge_is_subsgm(trg->s1))
					if (!vcn_container_exist(encroached_trg, trg->t1))
						if (!vcn_container_exist(unencroached_trg, trg->t1))
							vcn_container_insert(processing_trg, trg->t1);
			}
			if (trg->t2 != NULL) {
				if (!medge_is_subsgm(trg->s2))
					if (!vcn_container_exist(encroached_trg, trg->t2))
						if (!vcn_container_exist(unencroached_trg, trg->t2))
							vcn_container_insert(processing_trg, trg->t2);
			}
			if (trg->t3 != NULL){
				if (!medge_is_subsgm(trg->s3))
					if (!vcn_container_exist(encroached_trg, trg->t3))
						if (!vcn_container_exist(unencroached_trg, trg->t3))
							vcn_container_insert(processing_trg, trg->t3);
			}
		} else {
			vcn_container_insert(unencroached_trg, trg);
		}
	}
	/* Free memory */
	vcn_container_destroy(unencroached_trg);
	vcn_container_destroy(processing_trg);

	return encroached_trg;
}

static inline vcn_container_t* remove_encroached_triangles
             (vcn_mesh_t *const restrict mesh,
	      msh_trg_t *const restrict first_trg_to_check,
	      const msh_vtx_t *const restrict v,
	      /* v_orfans_reference: NULL to take an arbitrary reference */
	      const msh_vtx_t *const restrict v_orfans_reference,
	      /* big_trg: NULL if not required */
	      vcn_container_t *const restrict big_trg,
	      /* poor_quality_trg: NULL if not required */
	      hash_trg_t *const restrict poor_quality_trg)
{
	vcn_container_t *const restrict encroached_trg =
		get_encroached_triangles(first_trg_to_check, v);

	/* Destroy encroached triangles */
	vcn_container_t *orfan_vtx = vcn_container_create(VCN_CONTAINER_SORTED);
	while (vcn_container_is_not_empty(encroached_trg)) {
		msh_trg_t *const restrict trg =
			vcn_container_delete_first(encroached_trg);
		/* Set vertices as orfans :( */
		vcn_container_insert(orfan_vtx, trg->v1);
		vcn_container_insert(orfan_vtx, trg->v2);
		vcn_container_insert(orfan_vtx, trg->v3);

		mesh_substract_triangle(mesh, trg);

		/* Substract from poor-quality set (if exist) */
		if (NULL != poor_quality_trg)
			hash_trg_remove(poor_quality_trg, trg);

		/* Substract from big triangles set (if exist) */
		if (NULL != big_trg)
			remove_big_trg(big_trg, trg);

		if (NULL != trg->attr)
			free(trg->attr);
		free(trg);
	}
	vcn_container_destroy(encroached_trg);

	/* Return orfan vertices */
	return orfan_vtx;
}

static inline msh_vtx_t* retrg_fan_get_next_trg
                          (uint32_t* N_orfan_vtx,
			   msh_vtx_t ** orfan_vtx,
			   const msh_vtx_t *const v1,
			   const msh_vtx_t *const v2){
  double reference_angle = atan2(v2->x[1] - v1->x[1],
				 v2->x[0] - v1->x[0]);
  double min_angle = VCN_MATH_PI;
  uint32_t i3 = *N_orfan_vtx;
  for (uint32_t i = 0; i < *N_orfan_vtx; i++) {
    if (v2 == orfan_vtx[i])
      continue;
    double angle = atan2(orfan_vtx[i]->x[1] - v1->x[1],
			 orfan_vtx[i]->x[0] - v1->x[0]);
    double angle_diff = angle - reference_angle;
    if (angle_diff < 0.0)
      angle_diff += 2.0 * VCN_MATH_PI;
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

static inline void retriangulate_fan
                             (vcn_mesh_t *const restrict mesh,
			      const msh_vtx_t *const restrict v_pivot,
			      const msh_vtx_t *const restrict v_start,
			      vcn_container_t *restrict orfan_vtx,
			      /* l_new_trg: NULL if not required */
			      vcn_container_t *const restrict l_new_trg)
{
	/* Allocate orfan vertices into an array */
	msh_vtx_t** orfan_array = malloc(vcn_container_get_length(orfan_vtx)
				       * sizeof(*orfan_array));
	uint32_t N_orfan_vtx = 0;
	while (vcn_container_is_not_empty(orfan_vtx)) {
		msh_vtx_t *restrict vtx = vcn_container_delete_first(orfan_vtx);
		orfan_array[N_orfan_vtx] = vtx;
		N_orfan_vtx += 1;
	}

	/* Retriangulate fan */
	msh_vtx_t *restrict vtx = (msh_vtx_t*) v_start;
	msh_edge_t *restrict sgm;
	do {
		msh_vtx_t *restrict v3 = 
			retrg_fan_get_next_trg(&N_orfan_vtx, 
					       orfan_array,
					       v_pivot, vtx);

		/* Create triangle */
		msh_trg_t *restrict trg = mtrg_create();
		trg->v1 = (msh_vtx_t*) v_pivot;
		trg->v2 = vtx;
		trg->v3 = v3;

		/* Add triangle to the mesh */
		mesh_add_triangle(mesh, trg);
		mesh->do_after_insert_trg(mesh);

		if (l_new_trg != NULL)
			vcn_container_insert(l_new_trg, trg);

		vtx = v3;

		sgm = mesh_exist_edge(mesh->ht_edge, v_pivot, vtx);
	} while (vtx != v_start && !medge_is_subsgm(sgm));
	free(orfan_array);
}

static inline void verify_new_encroachments
                             (vcn_mesh_t *const restrict mesh,
			      const msh_vtx_t *const restrict v,
			      vcn_container_t *const restrict l_new_trg,
			      /* Could be NULL if not required */
			      vcn_container_t *const restrict encroached_sgm,
			      vcn_container_t *const restrict big_trg,
			      hash_trg_t *const restrict poor_quality_trg)
{
	bool (*inside)(const double[2], const double[2], const double [2]) =
		vcn_utils2D_pnt_lies_strictly_in_diametral_circle;

	while (vcn_container_is_not_empty(l_new_trg)) {
		msh_trg_t* restrict trg = vcn_container_delete_first(l_new_trg);
		if (NULL != encroached_sgm) {
			msh_edge_t* restrict edge =
				mtrg_get_opposite_edge(trg, v);
			if (medge_is_subsgm(edge)) {
				if (inside(edge->v1->x, edge->v2->x, v->x)) {
					vcn_container_insert(encroached_sgm,
							     edge);
					continue;
				}
			}
		}		
		check_trg(trg, mesh, big_trg, poor_quality_trg);
	}
}

void insert_vertex(vcn_mesh_t *const restrict mesh,
		   msh_trg_t *const restrict trg_containing_cc,
		   const msh_vtx_t *const restrict cc,
		   /* big_trg: NULL if not required */
		   vcn_container_t *const restrict  big_trg,
		   /* poor_quality_trg: NULL if not required */
		   hash_trg_t *const restrict poor_quality_trg,
		   /* l_new_trg: NULL if not required */
		   vcn_container_t *const restrict l_new_trg)
{
	/* Insert vertex in the built-in structure */
	vcn_bins2D_insert(mesh->ug_vtx, cc);

	/* Remove encroached triangles */
	vcn_container_t* orfan_vertices =
		remove_encroached_triangles(mesh, trg_containing_cc, 
					    cc, NULL, big_trg,
					    poor_quality_trg);
	retriangulate_fan(mesh, cc, 
			  vcn_container_get_first(orfan_vertices),
			  orfan_vertices,
			  l_new_trg);
	vcn_container_destroy(orfan_vertices);
	mesh->do_after_insert_vtx(mesh);
}

static inline void insert_midpoint
                             (vcn_mesh_t *const restrict mesh,
			      msh_edge_t *const restrict sgm,
			      const msh_vtx_t *const restrict v,
			      vcn_container_t *const restrict big_trg,
			      hash_trg_t *const restrict poor_quality_trg,
			      msh_edge_t* subsgm[2],
			      vcn_container_t *const restrict l_new_trg)
{
	/* Insert vertex in the built-in structure */
	vcn_bins2D_insert(mesh->ug_vtx, v);

	/* Split the subsegment */  
	subsgm[0] = (msh_edge_t*) calloc(1, sizeof(msh_edge_t));
	subsgm[0]->v1 = sgm->v1;
	subsgm[0]->v2 = (msh_vtx_t*) v;
	vcn_container_insert(mesh->ht_edge, subsgm[0]);
	
	subsgm[1] = (msh_edge_t*) calloc(1, sizeof(msh_edge_t));
	subsgm[1]->v1 = (msh_vtx_t*) v;
	subsgm[1]->v2 = sgm->v2;
	vcn_container_insert(mesh->ht_edge, subsgm[1]);

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
		vcn_container_t* orfan_vertices =
			remove_encroached_triangles(mesh, sgm->t1, v, sgm->v2,
						    big_trg, poor_quality_trg);
		retriangulate_fan(mesh, v, sgm->v2, orfan_vertices, l_new_trg);
		vcn_container_destroy(orfan_vertices);
	}
	/* Retriangulate right */
	if (NULL != sgm->t2) {
		vcn_container_t* orfan_vertices =
			remove_encroached_triangles(mesh, sgm->t2,
						    v, sgm->v1,
						    big_trg,
						    poor_quality_trg);
		retriangulate_fan(mesh, v, sgm->v1, orfan_vertices, l_new_trg);
		vcn_container_destroy(orfan_vertices);
	}
	/* Free splitted segment */
	vcn_container_delete(mesh->ht_edge, sgm);
	medge_destroy_subsgm_attribute(sgm);
	free(sgm);
}

static inline void split_encroached_segments
                             (vcn_mesh_t *const restrict mesh,
			      vcn_container_t *const restrict encroached_sgm,
			      vcn_container_t *const restrict big_trg,
			      hash_trg_t *const restrict poor_quality_trg)
{
	while (vcn_container_is_not_empty(encroached_sgm)) {
		msh_edge_t* sgm = vcn_container_delete_first(encroached_sgm);
		msh_vtx_t *v = get_midpoint(sgm);
		/* Split the segment */
		vcn_container_t* l_new_trg = vcn_container_create(VCN_CONTAINER_QUEUE);
		msh_edge_t* subsgm[2];
		insert_midpoint(mesh, sgm, v, big_trg, poor_quality_trg,
				subsgm, l_new_trg);

		/* Verify new encroachments */
		verify_new_encroachments(mesh, v, l_new_trg, encroached_sgm,
					 big_trg, poor_quality_trg);
		vcn_container_destroy(l_new_trg);

		/* Check if the first subsegment is encroached */
		if (is_encroached(subsgm[0]))
			vcn_container_insert(encroached_sgm, subsgm[0]);

		/* Check if the second subsegment is encroached */
		if (is_encroached(subsgm[1]))
			vcn_container_insert(encroached_sgm, subsgm[1]);

		mesh->do_after_insert_vtx(mesh);
	}
}

static inline void initialize_encroached_sgm
                          (vcn_mesh_t *const restrict mesh,
			   vcn_container_t *const restrict encroached_sgm)
{
	for (uint32_t i = 0; i < mesh->N_input_sgm; i++) {
		msh_edge_t* restrict sgm = mesh->input_sgm[i];
		while (sgm != NULL) {
			if (is_encroached(sgm))
				vcn_container_insert(encroached_sgm, sgm);
			sgm = medge_subsgm_next(sgm);
		}
	}
}

static inline void insert_big_trg
                         (vcn_container_t * restrict big_trg,
			  msh_trg_t * restrict trg,
			  double big_ratio)
{
	uint32_t *attr = malloc(sizeof(uint32_t));
	*attr = (uint32_t) (1e2 / big_ratio);
	trg->attr = attr;
	vcn_container_insert(big_trg, trg);
}

static inline msh_trg_t* remove_bigger_trg(vcn_container_t * restrict big_trg)
{
	msh_trg_t *trg = vcn_container_delete_first(big_trg);
	if (NULL != trg) {
		free(trg->attr);
		trg->attr = NULL;
	}
	return trg;  
}

static inline void remove_big_trg(vcn_container_t *big_trg, msh_trg_t *trg)
{
	if (NULL != trg->attr) {
		if (trg == vcn_container_delete(big_trg, trg)) {
			free(trg->attr);
			trg->attr = NULL;
		}
	}
}

static inline void check_trg(msh_trg_t *const restrict trg,
			     vcn_mesh_t *const restrict mesh,
			     vcn_container_t *const restrict big_trg,
			     hash_trg_t *const restrict poor_quality_trg)
{
	double big_ratio;
	if (mtrg_is_too_big(mesh, trg, &big_ratio)) {
		insert_big_trg(big_trg, trg, big_ratio);
	} else {
		if (1e30 > mesh->cr2se_ratio) {
			double trg_cr2se_ratio =
				vcn_utils2D_get_cr2se_ratio(trg->v1->x,
							    trg->v2->x,
							    trg->v3->x);
			if (trg_cr2se_ratio > mesh->cr2se_ratio)
				hash_trg_insert(poor_quality_trg,
						trg, trg_cr2se_ratio);
		}
	}
}

static void initialize_big_and_poor_quality_trg
                          (vcn_mesh_t *const restrict mesh,
			   vcn_container_t *const restrict big_trg,
			   hash_trg_t *const restrict poor_quality_trg)
{
	vcn_iterator_t *const restrict iter = vcn_iterator_create();
	vcn_iterator_set_container(iter, mesh->ht_trg);
	while (vcn_iterator_has_more(iter)) {
		msh_trg_t *restrict trg =
			(msh_trg_t*) vcn_iterator_get_next(iter);
		check_trg(trg, mesh, big_trg, poor_quality_trg);
	}
	vcn_iterator_destroy(iter);
}

static inline msh_trg_t* get_trg_containing_circumcenter
                          (const vcn_mesh_t *const restrict mesh,
			   const msh_trg_t *const restrict trg,
			   const msh_vtx_t *const restrict cc)
{
	/* Walk direct to the circumcenter */
	msh_trg_t *restrict prev2_trg = NULL;
	msh_trg_t *restrict prev1_trg = NULL;
	msh_trg_t *restrict trg_containing_cc = (msh_trg_t*) trg;
	while (!vcn_utils2D_pnt_lies_in_trg(trg_containing_cc->v1->x,
						    trg_containing_cc->v2->x,
						    trg_containing_cc->v3->x,
						    cc->x)) {
		/* Update history */
		prev2_trg = prev1_trg;
		prev1_trg = trg_containing_cc;

		/* Calculate centroid of the triangle */
		double centroid[2];
		vcn_utils2D_get_trg_centroid(trg_containing_cc->v1->x,
					     trg_containing_cc->v2->x,
					     trg_containing_cc->v3->x,
					     centroid);

		/* Move to the next triangle */
		if (vcn_utils2D_are_sgm_intersected(centroid, cc->x,
						    trg_containing_cc->v1->x,
						    trg_containing_cc->v2->x,
						    NULL, NULL))
			trg_containing_cc = trg_containing_cc->t1;
		else if (vcn_utils2D_are_sgm_intersected(centroid, cc->x,
							 trg_containing_cc->v2->x,
							 trg_containing_cc->v3->x,
							 NULL, NULL))
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

static inline vcn_container_t* get_sgm_encroached_by_vertex
                          (const vcn_mesh_t *const restrict mesh,
			   const msh_trg_t *const restrict trg_containing_vtx,
			   const msh_vtx_t *const restrict vtx)
{
	/* Get encroached triangles */
	vcn_container_t *const restrict encroached_trg =
		get_encroached_triangles(trg_containing_vtx, vtx);

	/* Get segments */
	vcn_container_t *const restrict segments =
		vcn_container_create(VCN_CONTAINER_SORTED);
	while (vcn_container_is_not_empty(encroached_trg)) {
		msh_trg_t *const restrict trg = 
			vcn_container_delete_first(encroached_trg);
		if (medge_is_subsgm(trg->s1))
			vcn_container_insert(segments, trg->s1);
		if (medge_is_subsgm(trg->s2))
			vcn_container_insert(segments, trg->s2);
		if (medge_is_subsgm(trg->s3))
			vcn_container_insert(segments, trg->s3);
	}
	/* Free memory */
	vcn_container_destroy(encroached_trg);

	/* Get encroached segments */
	vcn_container_t *restrict encroached_sgm =
		vcn_container_create(VCN_CONTAINER_QUEUE);
	while (vcn_container_is_not_empty(segments)) {
		msh_edge_t *const restrict sgm = 
			vcn_container_delete_first(segments);
		if (vcn_utils2D_pnt_lies_strictly_in_diametral_circle(sgm->v1->x,
								      sgm->v2->x,
								      vtx->x))
			vcn_container_insert(encroached_sgm, sgm);
	}

	/* Free memory */
	vcn_container_destroy(segments);

	return encroached_sgm;
}

static inline bool split_is_permitted(const msh_edge_t *const restrict sgm,
				      double d)
{
	/* Check if the length is not a power of two */
	double length = vcn_utils2D_get_dist(sgm->v1->x, sgm->v2->x);
	double k_real = vcn_math_log2(length / _VCN_CC_SHELL_UNIT);
	int k_int = (int)((k_real > 0)?(k_real + 0.5):(k_real - 0.5));
	double shell_length = _VCN_CC_SHELL_UNIT * pow(2.0, k_int);
	if (fabs(length - shell_length) > VCN_GEOMETRIC_TOL)
		return true;

	/* Get subsegment clusters */
	double smallest_angle;
	vcn_container_t* cluster1 =
		get_subsgm_cluster(sgm, sgm->v1, &smallest_angle);

	double smallest_angle_aux;
	vcn_container_t* cluster2 =
		get_subsgm_cluster(sgm, sgm->v2, &smallest_angle_aux);

	/* Check if s belongs to no cluster or if sgm belongs to 
	 * both clusters */
	if (vcn_container_is_empty(cluster1) == 
	    vcn_container_is_empty(cluster2)) {
		vcn_container_destroy(cluster1);
		vcn_container_destroy(cluster2);
		return true;
	}

	/* Get the cluster containing the sgm */
	vcn_container_t *restrict cluster;
	if (vcn_container_is_not_empty(cluster1)) {
		cluster = cluster1;
		vcn_container_destroy(cluster2);
	} else {
		cluster = cluster2;
		vcn_container_destroy(cluster1);
		smallest_angle = smallest_angle_aux;
	}

	/* Check if the cluster contains a segment shorter than sgm */
	while (vcn_container_is_not_empty(cluster)) {
		const msh_edge_t *const restrict cluster_sgm = 
			vcn_container_delete_first(cluster);

		if (cluster_sgm == sgm) 
			continue;
		double length_aux = vcn_utils2D_get_dist(cluster_sgm->v1->x,
					     cluster_sgm->v2->x);
		if (length - length_aux > VCN_GEOMETRIC_TOL) {
			vcn_container_destroy(cluster);
			return true;
		}
	}
	vcn_container_destroy(cluster);
  
	/* No new edge will be shorter than shortest existing edge */
	double r_min = length * sin(smallest_angle / 2.0);
	return r_min >= d + VCN_GEOMETRIC_TOL;
}

static inline vcn_container_t* get_subsgm_cluster 
                          (const msh_edge_t *const sgm,
			   const msh_vtx_t *const vtx,
			   double* smallest_angle)
{
	/* Create cluster */
	vcn_container_t* cluster = vcn_container_create(VCN_CONTAINER_SORTED);
	*smallest_angle = VCN_MATH_PI;

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
			angle += 2.0 * VCN_MATH_PI;
    
		/* If the angle is lower than 60 degrees, add to cluster */
		if (angle < VCN_MATH_PI / 3.0) {
			vcn_container_insert(cluster, sgm_next);
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
			angle += 2.0 * VCN_MATH_PI;
    
		/* If the angle is lower than 60 degrees, add to cluster */
		if (angle < VCN_MATH_PI / 3.0) {
			vcn_container_insert(cluster, sgm_next);
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
	if (vcn_container_is_not_empty(cluster))
		vcn_container_insert(cluster, sgm);
	else
		*smallest_angle = 0.0;
	return cluster;
}

static inline bool has_edge_length_constrained(const vcn_mesh_t *const mesh)
{
	return mesh->max_edge_length * mesh->scale > VCN_GEOMETRIC_TOL ||
		mesh->max_subsgm_length * mesh->scale > VCN_GEOMETRIC_TOL;
}

static bool edge_violates_constrain(const vcn_mesh_t *const restrict mesh,
				    const msh_edge_t *const restrict sgm,
				    /* big_ratio can be NULL if not required */
				    double *big_ratio)
{
	bool violates_constrain = false;
	double d = vcn_utils2D_get_dist(sgm->v1->x, sgm->v2->x);
	double edge_max = mesh->max_edge_length * mesh->scale;
	if (edge_max > VCN_GEOMETRIC_TOL && edge_max < d) {
		if (NULL != big_ratio)
			*big_ratio = d;
		violates_constrain = true;
	} else if (medge_is_subsgm(sgm)) {
		double sgm_max = mesh->max_subsgm_length * mesh->scale;
		if (sgm_max > VCN_GEOMETRIC_TOL && sgm_max < d) {
			if (NULL != big_ratio)
				*big_ratio = d;
			violates_constrain = true;
		}
	}
	return violates_constrain;
}

static bool edge_greater_than_density(const vcn_mesh_t *const restrict mesh,
				      const msh_edge_t *const restrict sgm,
				      /* big_ratio can be NULL if not required */
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
	double d = vcn_utils2D_get_dist(sgm->v1->x, sgm->v2->x);
	if (VCN_GEOMETRIC_TOL > d) {
		/* Prevent of eternal running */
		if (NULL != big_ratio)
			*big_ratio = VCN_GEOMETRIC_TOL;
		is_greater = false;  
	} else {
		d /= mesh->scale;
		double lh = calculate_lh(mesh, sgm, d, 1);
		if (_VCN_MAX_LH_TOLERATED < lh) {
			is_greater = true;
		} else {
			lh = calculate_lh(mesh, sgm, d, 2);
			is_greater = (_VCN_MAX_LH_TOLERATED < lh);
		}
		if (NULL != big_ratio)
			*big_ratio = lh;
	}
	return is_greater;
}

static double calculate_lh(const vcn_mesh_t *const restrict mesh,
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
		if (h1 > _VCN_MAX_GRADING_RATIO * h2)
			h1 = _VCN_MAX_GRADING_RATIO * h2;
		else if (h2 > _VCN_MAX_GRADING_RATIO * h1)
			h2 = _VCN_MAX_GRADING_RATIO * h1;
		
		integral += (h1 + h2);
		h1 = h2;
	}
	return 2 * dist / (integral * width);
}

static inline double calculate_h(const vcn_mesh_t *const restrict mesh,
				 const msh_edge_t *const restrict sgm,
				 double t)
{
	double vtx[2];
	vtx[0] = (1 - t) * sgm->v1->x[0] + t * sgm->v2->x[0];
	vtx[1] = (1 - t) * sgm->v1->x[1] + t * sgm->v2->x[1];
	double vtx_scaled[2];
	mesh_get_extern_scale_and_disp(mesh, vtx, vtx_scaled);
	double density = mesh->density(vtx_scaled, mesh->density_data);
	if (density < VCN_GEOMETRIC_TOL)
		density = VCN_GEOMETRIC_TOL;
	return 1.0 / density;
}

static inline bool mtrg_is_too_big(const vcn_mesh_t *const restrict mesh,
				   const msh_trg_t *const restrict trg,
				   /* big_ratio could be NULL if not required */
				   double *big_ratio)
{
	if(NULL != big_ratio)
		*big_ratio = 1.0;
	bool is_too_big = false;
	if (has_edge_length_constrained(mesh)) {
		is_too_big = edge_violates_constrain(mesh, trg->s1, big_ratio);
		if (!is_too_big)
			is_too_big = edge_violates_constrain(mesh, trg->s2,
							     big_ratio);
		if (!is_too_big)
			is_too_big = edge_violates_constrain(mesh, trg->s3,
							     big_ratio);
	} else {
		if (NULL != mesh->density) {
			is_too_big = edge_greater_than_density(mesh, trg->s1,
							       big_ratio);
			if (!is_too_big)
				is_too_big = 
					edge_greater_than_density(mesh,
								  trg->s2,
								  big_ratio);
			if (!is_too_big)
				is_too_big = 
					edge_greater_than_density(mesh,
								  trg->s3,
								  big_ratio);
		}
	}
	return is_too_big;
}
