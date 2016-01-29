/******************************************************************************
 *   Geometric Bot: Geometric tesselations for numerical analysis.            *
 *   2011-2015 Victor Eduardo Cardoso Nungaray                                *
 *   Twitter: @victore_cardoso                                                *
 *   email: victorc@cimat.mx                                                  *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "vcn/math_bot.h"
#include "vcn/container_bot.h"
#include "vcn/eigen_bot.h"
#include "vcn/geometric_bot/point2D.h"
#include "vcn/geometric_bot/utils2D.h"
#include "vcn/geometric_bot/knn/bins2D.h"
#include "vcn/geometric_bot/knn/bins2D_iterator.h"
#include "vcn/geometric_bot/model/model2D.h"
#include "vcn/geometric_bot/mesh/constrained_delaunay.h"
#include "vcn/geometric_bot/mesh/density2D.h"
#include "vcn/geometric_bot/mesh/mesh2D.h"

#include "../model/model2D_struct.h"
#include "mesh2D_structs.h"

#define _VCN_MAX_LH_TOLERATED (1.5)
#define _VCN_MAX_GRADING_RATIO (27.0)
#define _VCN_SUBSEGMENT_VTX ((void*)0x2)
#define _VCN_CC_SHELL_UNIT (1e-3)
   
/******************** Private structures *****************************/
typedef struct {
	vcn_container_t* avl[64];
	uint32_t length;
} hash_trg_t;

/**************** Private functions **********************/
static hash_trg_t* hash_trg_create(void);

static uint32_t hash_trg_length(const hash_trg_t *const htrg);

static void hash_trg_insert(hash_trg_t *const htrg,
				   msh_trg_t *const trg,
				   double cr2se_ratio);
static msh_trg_t* hash_trg_remove_first(hash_trg_t *const htrg);

static void hash_trg_remove(hash_trg_t *const htrg,
				   msh_trg_t *const trg);

static void hash_trg_destroy(hash_trg_t *const htrg);

static bool medge_is_too_big(const vcn_mesh_t *const mesh,
			     const msh_edge_t *const sgm,
			     double (*density)(const void *const data,
					       const double const x[2]),
			     const void *const density_data,
			     /* big_ratio can be NULL if not required */
			     double *big_ratio);

static bool mtrg_is_too_big(const vcn_mesh_t *const restrict mesh,
			    const msh_trg_t *const restrict trg,
			     double (*density)(const void *const data,
					       const double const x[2]),
			    const void *const restrict density_data,
			    /* big_ratio could be NULL if not required */
			    double *big_ratio);

static void mesh_get_extern_scale_and_disp
                           (const vcn_mesh_t *const mesh,
			    const double *const internal,
			    double external[2]);

static void mesh_reallocate_htables(vcn_mesh_t *const mesh);

static msh_trg_t* mesh_locate_vtx(const vcn_mesh_t *const mesh,
				   const msh_vtx_t *const v);

static void null_do_after_trg(const vcn_mesh_t *const mesh);

/* Functions to handle the hash table storing the segments */
static bool refinement_check_max_vtx(const vcn_mesh_t *const mesh,
					    uint32_t max_vtx);

static bool refinement_check_max_trg(const vcn_mesh_t *const mesh,
					    uint32_t max_trg);

static bool refinement_is_encroached(const msh_edge_t *const sgm);

static msh_vtx_t* refinement_get_midpoint
                         (const msh_edge_t *const sgm);

static vcn_container_t* refinement_get_encroached_triangles
                             (const msh_trg_t *const first_trg_to_check,
			      const msh_vtx_t *const v);

static vcn_container_t* refinement_remove_encroached_triangles
             (vcn_mesh_t *const mesh,
	      msh_trg_t *const first_trg_to_check,
	      const msh_vtx_t *const v,
	      /* v_orfans_reference: NULL to take an arbitrary reference */
	      const msh_vtx_t *const v_orfans_reference,
	      /* big_trg: NULL if not required */
	      vcn_container_t *const big_trg,
	      /* poor_quality_trg: NULL if not required */
	      hash_trg_t *const poor_quality_trg);

static msh_vtx_t* refinement_retrg_fan_get_next_trg
                          (uint32_t* N_orfan_vtx,
			   msh_vtx_t ** orfan_vtx,
			   const msh_vtx_t *const v1,
			   const msh_vtx_t *const v2);

static void refinement_retriangulate_fan
                             (vcn_mesh_t *const mesh,
			      const msh_vtx_t *const v_pivot,
			      const msh_vtx_t *const v_start,
			      vcn_container_t *orfan_vtx,
			      /* l_new_trg: NULL if not required */
			      vcn_container_t *const l_new_trg);

static void refinement_insert_midpoint
                             (vcn_mesh_t *const mesh,
			      msh_edge_t *const sgm,
			      const msh_vtx_t *const v,
			      vcn_container_t *const big_trg,
			      hash_trg_t *const poor_quality_trg,
			      msh_edge_t* subsgm[2],
			      vcn_container_t *const l_new_trg);

static void refinement_insert_circumcenter
                             (vcn_mesh_t *const mesh,
			      msh_trg_t *const trg_containing_cc,
			      const msh_vtx_t *const cc,
			      /* big_trg: NULL if not required */
			      vcn_container_t *const big_trg,
			      /* poor_quality_trg: NULL if not required */
			      hash_trg_t *const poor_quality_trg,
			      /* l_new_trg: NULL if not required */
			      vcn_container_t *const l_new_trg);

static void refinement_split_encroached_segments
                                   (vcn_mesh_t *const mesh,
				    vcn_container_t *const encroached_sgm,
				    vcn_container_t *const big_trg,
				    hash_trg_t *const poor_quality_trg,
				    double cr2se_ratio,
				    double (*density)(const void *const data,
						      const double const x[2]),
				    const void *const density_data);
						   
static void refinement_verify_new_encroachments
                             (vcn_mesh_t *const mesh,
			      const msh_vtx_t *const  v,
			      vcn_container_t *const l_new_trg,
			      vcn_container_t *const encroached_sgm,
			      vcn_container_t *const big_trg,
			      hash_trg_t *const poor_quality_trg,
			      double cr2se_ratio,
			      double (*density)(const void *const data,
						const double const x[2]),
			      const void *const density_data);

static void refinement_initialize_encroached_sgm
                          (vcn_mesh_t *const mesh,
			   vcn_container_t *const encroached_sgm);

static void refinement_insert_big_trg
                   (vcn_container_t *big_trg, msh_trg_t *trg, double big_ratio);

static msh_trg_t* refinement_remove_bigger_trg(vcn_container_t *big_trg);

static void refinement_remove_big_trg
                             (vcn_container_t *big_trg, msh_trg_t *trg);

static void refinement_check_trg
                          (msh_trg_t *const trg, vcn_mesh_t *const mesh,
			   vcn_container_t *const big_trg,
			   hash_trg_t *const poor_quality_trg,
			   double cr2se_ratio,
			   double (*density)(const void *const data,
					     const double const x[2]),
			   const void *const density_data);

static void refinement_initialize_big_and_poor_quality_trg
                          (vcn_mesh_t *const mesh,
			   vcn_container_t *const big_trg,
			   hash_trg_t *const poor_quality_trg,
			   double cr2se_ratio,
			   double (*density)(const void *const data,
					     const double const x[2]),
			   const void *const density_data);

static msh_trg_t* refinement_get_trg_containing_circumcenter
                          (const vcn_mesh_t *const mesh,
			   const msh_trg_t *const trg,
			   const msh_vtx_t *const cc);

static vcn_container_t* refinement_get_sgm_encroached_by_vertex
                          (const vcn_mesh_t *const mesh,
			   const msh_trg_t *const trg_containing_vtx,
			   const msh_vtx_t *const vtx);

static bool refinement_split_is_permitted
                          (const msh_edge_t *const sgm,
			   double d);

static vcn_container_t* refinement_get_subsgm_cluster 
                          (const msh_edge_t *const sgm,
			   const msh_vtx_t *const vtx,
			   double* smallest_angle);

/* Triangulation functions */
static void remove_triangles_propagate
                          (vcn_mesh_t *const mesh,
			   msh_trg_t* trg);

static void remove_concavities_triangles(vcn_mesh_t* mesh);
static void remove_holes_triangles(vcn_mesh_t* mesh, double* holes, uint32_t N_holes);

/* Compare functions */
static int compare_trg_attr_uint64_t(const void *const trgA,
				     const void *const trgB);

static int compare_deterministic_sgm(const void *const sgmA,
				     const void *const sgmB);
static int compare_sgmA_isSmallerThan_sgmB(const void *const sgmA, 
					   const void *const sgmB);
static int compare_sgmA_isTheSameThan_sgmB(const void *const  sgmA, 
					   const void *const  sgmB);
static int compare_trgA_isBetterThan_trgB(const void *const  trgA, 
					  const void *const  trgB);
static int compare_trgA_isSmallerThan_trgB(const void *const  trgA,
					   const void *const  trgB);
static int compare_area1_isGreaterThan_area2(const void *const  a1,
					     const void *const  a2);
static int compare_id(const void* const ptrA,
		      const void* const ptrB);
static int compare_vtx(const void* const vtxA,
		       const void* const vtxB);
static int compare_sgm_by_dist_and_by_vtx_pointer
				(const void* const sgmA,
				 const void* const sgmB);

/* Hash functions */
static uint32_t hash_key_vtx(const void *const  vertex);
static uint32_t hash_key_trg(const void *const  triangle);

static inline hash_trg_t* hash_trg_create(void)
{
	hash_trg_t* htrg = calloc(1, sizeof(hash_trg_t));
	for (uint32_t i = 0; i < 64; i++) {
		htrg->avl[i] = vcn_container_create(VCN_CONTAINER_SORTED);
		/* REFACTOR
		vcn_container_key_generator(htrg->avl[i], compare_trg_attr_uint64_t);
		vcn_container_comparer(htrg->avl[i], );
		*/
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
  if(trg->attr != NULL) return; /* It is already inserted */

  /* Calculate hash key */
  uint32_t hash_key = (uint32_t) cr2se_ratio;
  if (hash_key > 63)
    hash_key = 63;
  
  /* Set circumradius to shortest edge ratio as attribute */
  uint64_t* attr = malloc(sizeof(uint64_t));
  attr[0] = (uint64_t)(cr2se_ratio * 1e6);
  trg->attr = attr;

  /* Insert into the AVL */
  bool is_inserted =
    vcn_container_insert(htrg->avl[hash_key], trg);

  if (is_inserted) {
    htrg->length += 1;
  } else {
    free(trg->attr);
    trg->attr = NULL;
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
  if(trg->attr == NULL) return; /* It is not inserted */

  /* Calculate hash key */
  double cr2se_ratio = ((uint64_t*)trg->attr)[0] / 1e6;
  uint32_t hash_key = (uint32_t) cr2se_ratio;
  if(hash_key > 63) hash_key = 63;
  
  /* Remove from AVL if the trg exist */
  bool is_removed = vcn_container_delete(htrg->avl[hash_key], trg);

  if (is_removed) {
    htrg->length -= 1;
    free(trg->attr);
    trg->attr = NULL;
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

static bool medge_is_too_big(const vcn_mesh_t *const restrict mesh,
			     const msh_edge_t *const restrict sgm,
			     double (*density)(const void *const data,
					       const double const x[2]),
			     const void *const restrict density_data,
			     /* big_ratio can be NULL if not required */
			     double *big_ratio)
{
	/* Get length of the segment (scaled) */
	double d = vcn_utils2D_get_dist(sgm->v1->x, sgm->v2->x);

	/* REFACTOR next line */
	if (density == VCN_DENSITY_MAX) {
		/* Constraint max values */
		double* restrict data = (double*) density_data;
		double edge_max = data[0] * mesh->scale;
		if (edge_max > VCN_GEOMETRIC_TOL && edge_max < d) {
			if (NULL != big_ratio)
				*big_ratio = d;
			return true;
		}

		if (medge_is_subsgm(sgm)) {
			double sgm_max = data[1] * mesh->scale;
			if (sgm_max > VCN_GEOMETRIC_TOL && sgm_max < d) {
				if (NULL != big_ratio)
					*big_ratio = d;
				return true;
			}
		}
		return false;
	}
	
	/* REFACTOR next line */
	if (density == VCN_DENSITY_IMG)
		/* The density is given by an image */
		density = vcn_density_img_get_density;

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
	 *
	 *    For convenience, we estimate
	 *
	 *                   1 _             x(t) = v1 + t * (v2-v1)
	 *              /     |              e.g.
	 *    lh =  1  /  |s| |  h(x) dt,          x( 0 ) = v1
	 *            /      _|                    x(0.5) = 0.5 * (v1 + v2)
	 *                   0                     x( 1 ) = v2
	 */

	/* Prevent of eternal running */
	if (d < VCN_GEOMETRIC_TOL) {
		if (NULL != big_ratio)
			*big_ratio = VCN_GEOMETRIC_TOL;
		return false;  
	}

	/* Scale distance */
	d /= mesh->scale;

	/* Integrate using one trapezoid */
	double v1_extern[2];
	mesh_get_extern_scale_and_disp(mesh, sgm->v1->x, v1_extern);
	double density_1 = density(v1_extern, density_data);
	if (density_1 < VCN_GEOMETRIC_TOL)
		density_1 = VCN_GEOMETRIC_TOL;
	double h1 = 1.0 / density_1;

	double v2_extern[2];
	mesh_get_extern_scale_and_disp(mesh, sgm->v2->x, v2_extern);
	double density_2 = density(v2_extern, density_data);
	if (density_2 < VCN_GEOMETRIC_TOL)
		density_2 = VCN_GEOMETRIC_TOL;
	double h2 = 1.0 / density_2;

	if (h1 > _VCN_MAX_GRADING_RATIO * h2)
		h1 = _VCN_MAX_GRADING_RATIO * h2;
	else if (h2 > _VCN_MAX_GRADING_RATIO * h1)
		h2 = _VCN_MAX_GRADING_RATIO * h1;

	double lh = d * 2.0 / (h1 + h2);

	if (lh > _VCN_MAX_LH_TOLERATED) {
		if (NULL != big_ratio)
			*big_ratio = lh;
		return true;
	}
  
	/* Integrate using two trapezoids */
	double midpoint[2];
	midpoint[0] = (sgm->v1->x[0] + sgm->v2->x[0]) / 2.0;
	midpoint[1] = (sgm->v1->x[1] + sgm->v2->x[1]) / 2.0;

	double mp_extern[2];
	mesh_get_extern_scale_and_disp(mesh, midpoint, mp_extern);
	double density_m = density(mp_extern, density_data);
	if (density_m < VCN_GEOMETRIC_TOL)
		density_m = VCN_GEOMETRIC_TOL;
	double h_midpoint = 1.0 / density_m;

	double min_h = h1;
	if (h2 < h1)
		min_h = h2;

	if (h_midpoint > _VCN_MAX_GRADING_RATIO * min_h) {
		h_midpoint = _VCN_MAX_GRADING_RATIO * min_h;
	} else if (h_midpoint < min_h) {
		if (h1 > _VCN_MAX_GRADING_RATIO * h_midpoint)
			h1 = _VCN_MAX_GRADING_RATIO * h_midpoint;
		if (h2 > _VCN_MAX_GRADING_RATIO * h_midpoint)
			h2 = _VCN_MAX_GRADING_RATIO * h_midpoint;
	}

	lh = d * 4.0 / (h1 + 2.0 * h_midpoint + h2);
	if (NULL != big_ratio)
		*big_ratio = lh;

	return lh > _VCN_MAX_LH_TOLERATED;
}

static void mesh_reallocate_htables(vcn_mesh_t *const restrict mesh)
{
	double avg_points_x_cell = 
		vcn_bins2D_get_length(mesh->ug_vtx) / 
		vcn_bins2D_get_N_bins(mesh->ug_vtx);
	double bin_size =
		vcn_bins2D_get_size_of_bins(mesh->ug_vtx) /
		sqrt(avg_points_x_cell);
  
	/* Allocate new grid */
	vcn_bins2D_t *ug_vtx = vcn_bins2D_create(bin_size);

	vcn_container_t *ht_edge = vcn_container_create(VCN_CONTAINER_HASH);
	vcn_container_set_key_generator(ht_edge, hash_key_edge);
	vcn_container_set_comparer(ht_edge, are_equal_edge);

	vcn_container_t *ht_trg = vcn_container_create(VCN_CONTAINER_HASH);
	vcn_container_set_key_generator(ht_trg, hash_key_trg);

	/* Reinsert vertices */
	while (vcn_bins2D_is_not_empty(mesh->ug_vtx)) {
		msh_vtx_t* vtx = vcn_bins2D_delete_first(mesh->ug_vtx);
		vcn_bins2D_insert(ug_vtx, vtx);
	}

	/* Point to new grid */
	vcn_bins2D_destroy(mesh->ug_vtx);
	mesh->ug_vtx = ug_vtx;

	/* Reinsert edges */
	while (vcn_container_is_not_empty(mesh->ht_edge)) {
		msh_edge_t* edge = vcn_container_delete_first(mesh->ht_edge);
		vcn_container_insert(ht_edge, edge);
	}
  
	/* Point to new hash table with edges */
	vcn_container_destroy(mesh->ht_edge);
	mesh->ht_edge = ht_edge;

	/* Reinsert triangles */
	while (vcn_container_is_not_empty(mesh->ht_trg)) {
		msh_trg_t* trg = vcn_container_delete_first(mesh->ht_trg);
		vcn_container_insert(ht_trg, trg);
	}
  
	/* Point to new hash table with triangles */
	vcn_container_destroy(mesh->ht_trg);
	mesh->ht_trg = ht_trg;
}

static inline bool mtrg_is_too_big(const vcn_mesh_t *const restrict mesh,
				   const msh_trg_t *const restrict trg,
				   double (*density)(const void *const data,
						     const double const x[2]),
				   const void *const restrict density_data,
				   /* big_ratio could be NULL if not required */
				   double *big_ratio)
{
  if(NULL != big_ratio)
    *big_ratio = 1.0;
  if (NULL == density) 
    return false;
  /* Check if at least one of its edges is too big */
  if (medge_is_too_big(mesh, trg->s1, density, density_data, big_ratio))
    return true;
  if (medge_is_too_big(mesh, trg->s2, density, density_data, big_ratio))
    return true;
  return medge_is_too_big(mesh, trg->s3, density, density_data, big_ratio);
}

static msh_trg_t* mesh_locate_vtx(const vcn_mesh_t *const restrict mesh,
				  const msh_vtx_t *const restrict v)
{
	vcn_iterator_t *iter = vcn_iterator_create();
	vcn_iterator_set_container(iter, mesh->ht_trg);
	msh_trg_t *enveloping_trg = NULL;
	while (vcn_iterator_has_more(iter)) {
		const msh_trg_t *trg = vcn_iterator_get_next(iter);
		if (vcn_utils2D_pnt_lies_in_trg(trg->v1->x,
							trg->v2->x,
							trg->v3->x,
							v->x)) {
			enveloping_trg = (msh_trg_t*) trg;
			break;
		}
	}
	vcn_iterator_destroy(iter);
	return enveloping_trg;
}

vcn_mesh_t* vcn_mesh_create(void)
{
	vcn_mesh_t *mesh = calloc(1, sizeof(*mesh));
	/* Initialize Mesh */
	mesh->ug_vtx = vcn_bins2D_create(1.0);

	mesh->ht_edge = vcn_container_create(VCN_CONTAINER_HASH);
	vcn_container_set_key_generator(mesh->ht_edge, hash_key_edge);
	vcn_container_set_comparer(mesh->ht_edge, are_equal_edge);

	mesh->ht_trg = vcn_container_create(VCN_CONTAINER_HASH);
	vcn_container_set_key_generator(mesh->ht_trg, hash_key_trg);
	mesh->scale = 1.0;

	mesh->do_after_trg = null_do_after_trg;
	return mesh;
}

static inline void null_do_after_trg(const vcn_mesh_t *const mesh)
{
	; /* Null statement */
}

void vcn_mesh_set_do_after_trg(vcn_mesh_t *mesh,
			       void (*do_after_trg)
				       	(const vcn_mesh_t *const))
{
	mesh->do_after_trg = do_after_trg;
}

void vcn_mesh_set_max_vtx(vcn_mesh_t* mesh, uint32_t N)
{
	; /* REFACTOR */
}

void vcn_mesh_set_max_trg(vcn_mesh_t* mesh, uint32_t N)
{
	; /* REFACTOR */
}

void vcn_mesh_set_min_angle(vcn_mesh_t* mesh, double angle)
{
	; /* REFACTOR */
}

void vcn_mesh_set_density_CDT(vcn_mesh_t* mesh)
{
	; /* REFACTOR */
}

void vcn_mesh_set_density_EDGE(vcn_mesh_t* mesh, double edge_max_size,
			       double subsgm_max_size)
{
	; /* REFACTOR */
}

void vcn_mesh_set_density_IMG(vcn_mesh_t* mesh, void *img_data)
{
	; /* REFACTOR */
}

void vcn_mesh_set_density_function(vcn_mesh_t* mesh,
				   double (*density)(double*))
{
	; /* REFACTOR */
}

bool vcn_mesh_is_empty(const vcn_mesh_t *const mesh)
{
	return true;/* REFACTOR */
}

static inline void mesh_get_extern_scale_and_disp
                           (const vcn_mesh_t *const mesh,
			    const double internal[2],
			    double external[2])
{
	external[0] = internal[0] / mesh->scale + mesh->xdisp;
	external[1] = internal[1] / mesh->scale + mesh->ydisp;
}



void vcn_mesh_destroy(vcn_mesh_t* mesh)
{
	free(mesh->input_vtx);
	if (mesh->N_input_sgm > 0) {
		/* Free attributes from the segments forming the input */
		for (uint32_t i = 0; i < mesh->N_input_sgm; i++) {
			msh_edge_t* sgm = mesh->input_sgm[i];
			while (NULL != sgm) {
				msh_edge_t* to_free = sgm;
				sgm = medge_subsgm_next(sgm);
				medge_destroy_subsgm_attribute(to_free);
			}
		}
		free(mesh->input_sgm);
	}
	/* Free mesh data */
	vcn_bins2D_destroy(mesh->ug_vtx);
	vcn_container_set_destroyer(mesh->ht_trg, free);
	vcn_container_destroy(mesh->ht_trg);
	vcn_container_set_destroyer(mesh->ht_edge, free);
	vcn_container_destroy(mesh->ht_edge);
	free(mesh);
}

bool vcn_mesh_is_vtx_inside(const vcn_mesh_t *const restrict mesh,
				 const double *const restrict vtx){
	msh_vtx_t v; 
	v.x[0] = mesh->scale * (vtx[0] - mesh->xdisp);
	v.x[1] = mesh->scale * (vtx[1] - mesh->ydisp);
	msh_trg_t* enveloping_trg = mesh_locate_vtx(mesh, &v);
	return (enveloping_trg != NULL);
}

void vcn_mesh_get_vertices(vcn_mesh_t* mesh, double* vertices)
{
	vcn_bins2D_iter_t* iter = vcn_bins2D_iter_create();
	vcn_bins2D_iter_set_bins(iter, mesh->ug_vtx);
	int i = 0;
	while (vcn_bins2D_iter_has_more(iter)) {
		const msh_vtx_t* vtx = vcn_bins2D_iter_get_next(iter);
		mesh_get_extern_scale_and_disp(mesh, vtx->x,
					       &(vertices[i * 2]));
		i++;
	}
	vcn_bins2D_iter_destroy(iter);
}

inline uint32_t vcn_mesh_get_N_vtx(const vcn_mesh_t *const mesh)
{
	return vcn_bins2D_get_length(mesh->ug_vtx);
}

inline uint32_t vcn_mesh_get_N_trg(const vcn_mesh_t *const mesh)
{
	return vcn_container_get_length(mesh->ht_trg);
}

inline uint32_t vcn_mesh_get_N_edg(const vcn_mesh_t *const mesh)
{
	return vcn_container_get_length(mesh->ht_edge);
}

double vcn_mesh_get_area(const vcn_mesh_t *const mesh)
{
	double area = 0.0;
	vcn_iterator_t* iter = vcn_iterator_create();
	vcn_iterator_set_container(iter, mesh->ht_trg);
	while (vcn_iterator_has_more(iter)) {
		msh_trg_t* trg = (msh_trg_t*) vcn_iterator_get_next(iter);
		area += vcn_utils2D_get_2x_trg_area(trg->v1->x,
						    trg->v2->x,
						    trg->v3->x);
	}
	vcn_iterator_destroy(iter);
	return (0.5 * area) / vcn_math_pow2(mesh->scale);
}

static inline bool refinement_check_max_vtx
                             (const vcn_mesh_t *const restrict mesh,
			      uint32_t max_vtx)
{
	return (max_vtx == 0)?true:(vcn_bins2D_get_length(mesh->ug_vtx) < max_vtx);
}

static inline bool refinement_check_max_trg
                             (const vcn_mesh_t *const restrict mesh,
			      uint32_t max_trg)
{
	return (max_trg == 0)?true:(vcn_container_get_length(mesh->ht_trg) < max_trg);
}

static inline bool refinement_is_encroached
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

static inline msh_vtx_t* refinement_get_midpoint
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

static inline vcn_container_t* refinement_get_encroached_triangles
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

static inline vcn_container_t* refinement_remove_encroached_triangles
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
	/* Get encroached triangles */
	vcn_container_t *const restrict encroached_trg =
		refinement_get_encroached_triangles(first_trg_to_check, v);

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
			refinement_remove_big_trg(big_trg, trg);

		free(trg);
	}
	vcn_container_destroy(encroached_trg);

	/* Return orfan vertices */
	return orfan_vtx;
}

static inline msh_vtx_t* refinement_retrg_fan_get_next_trg
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

static inline void refinement_retriangulate_fan
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
			refinement_retrg_fan_get_next_trg(&N_orfan_vtx, 
							  orfan_array,
							  v_pivot, vtx);

		/* Create triangle */
		msh_trg_t *restrict trg = mtrg_create();
		trg->v1 = (msh_vtx_t*) v_pivot;
		trg->v2 = vtx;
		trg->v3 = v3;

		/* Add triangle to the mesh */
		mesh_add_triangle(mesh, trg);

		if (l_new_trg != NULL)
			vcn_container_insert(l_new_trg, trg);

		vtx = v3;

		sgm = mesh_exist_edge(mesh->ht_edge, v_pivot, vtx);
	} while (vtx != v_start && !medge_is_subsgm(sgm));
	free(orfan_array);
}

static inline void refinement_verify_new_encroachments
                             (vcn_mesh_t *const restrict mesh,
			      const msh_vtx_t *const restrict v,
			      vcn_container_t *const restrict l_new_trg,
			      vcn_container_t *const restrict encroached_sgm,
			      vcn_container_t *const restrict big_trg,
			      hash_trg_t *const restrict poor_quality_trg,
			      double cr2se_ratio,
			      double (*density)(const void *const data,
						const double const x[2]),
			      const void *const restrict density_data)
{
  while (vcn_container_is_not_empty(l_new_trg)) {
    msh_trg_t* restrict trg = vcn_container_delete_first(l_new_trg);
    if (encroached_sgm != NULL) {
      msh_edge_t* restrict edge = mtrg_get_opposite_edge(trg, v);
      if (medge_is_subsgm(edge)) {
	if (vcn_utils2D_pnt_lies_strictly_in_diametral_circle(edge->v1->x,
							    edge->v2->x,
							    v->x)) {
	  vcn_container_insert(encroached_sgm, edge);
	  continue;
	}
      }
    }
    
    refinement_check_trg(trg, mesh, big_trg, poor_quality_trg,
			 cr2se_ratio, density, density_data);
  }
}

static inline void refinement_insert_midpoint
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
  if (prev != NULL)
    medge_update_subsgm_next(prev, subsgm[0]);
  if (next != NULL)
    medge_update_subsgm_prev(next, subsgm[1]);
  if (mesh->input_sgm[input_idx] == sgm)
    mesh->input_sgm[input_idx] = subsgm[0];
  
  /* Retriangulate left */
  if (sgm->t1 != NULL) {
    vcn_container_t* orfan_vertices =
      refinement_remove_encroached_triangles(mesh, sgm->t1, 
					     v, sgm->v2,
					     big_trg,
					     poor_quality_trg);
    refinement_retriangulate_fan(mesh, v, sgm->v2,
				 orfan_vertices,
				 l_new_trg);
    vcn_container_destroy(orfan_vertices);
  }
  /* Retriangulate right */
  if (sgm->t2 != NULL) {
    vcn_container_t* orfan_vertices =
      refinement_remove_encroached_triangles(mesh, sgm->t2,
					     v, sgm->v1,
					     big_trg,
					     poor_quality_trg);
    refinement_retriangulate_fan(mesh, v, sgm->v1,
				 orfan_vertices,
				 l_new_trg);
    vcn_container_destroy(orfan_vertices);
  }
  /* Free splitted segment */
  vcn_container_delete(mesh->ht_edge, sgm);
  medge_destroy_subsgm_attribute(sgm);
  free(sgm);
}

static inline void refinement_insert_circumcenter
                             (vcn_mesh_t *const restrict mesh,
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
    refinement_remove_encroached_triangles(mesh, trg_containing_cc, 
					   cc, NULL, big_trg,
					   poor_quality_trg);
  refinement_retriangulate_fan(mesh, cc, 
			       vcn_container_get_first(orfan_vertices),
			       orfan_vertices,
			       l_new_trg);
  vcn_container_destroy(orfan_vertices);
}

static inline void refinement_split_encroached_segments
                             (vcn_mesh_t *const restrict mesh,
			      vcn_container_t *const restrict encroached_sgm,
			      vcn_container_t *const restrict big_trg,
			      hash_trg_t *const restrict poor_quality_trg,
			      double cr2se_ratio,
			      double (*density)(const void *const data,
						const double const x[2]),
			      const void *const restrict density_data)
{
	while (vcn_container_is_not_empty(encroached_sgm)) {
		msh_edge_t* sgm = vcn_container_delete_first(encroached_sgm);
		msh_vtx_t *v = refinement_get_midpoint(sgm);
		/* Split the segment */
		vcn_container_t* l_new_trg = vcn_container_create(VCN_CONTAINER_QUEUE);
		msh_edge_t* subsgm[2];
		refinement_insert_midpoint(mesh, sgm, v, big_trg,
					   poor_quality_trg,
					   subsgm, l_new_trg);

		/* Verify new encroachments */
		refinement_verify_new_encroachments(mesh, v, l_new_trg, encroached_sgm,
						    big_trg, poor_quality_trg, cr2se_ratio,
						    density, density_data);
		vcn_container_destroy(l_new_trg);

		/* Check if the first subsegment is encroached */
		if (refinement_is_encroached(subsgm[0]))
			vcn_container_insert(encroached_sgm, subsgm[0]);

		/* Check if the second subsegment is encroached */
		if (refinement_is_encroached(subsgm[1]))
			vcn_container_insert(encroached_sgm, subsgm[1]);
	}
}

static inline void refinement_initialize_encroached_sgm
                          (vcn_mesh_t *const restrict mesh,
			   vcn_container_t *const restrict encroached_sgm)
{
	for (uint32_t i = 0; i < mesh->N_input_sgm; i++) {
		msh_edge_t* restrict sgm = mesh->input_sgm[i];
		while (sgm != NULL) {
			if(refinement_is_encroached(sgm))
				vcn_container_insert(encroached_sgm, sgm);
			sgm = medge_subsgm_next(sgm);
		}
	}
}

static inline void refinement_insert_big_trg
                         (vcn_container_t * restrict big_trg,
			  msh_trg_t * restrict trg,
			  double big_ratio)
{
	uint64_t *attr = malloc(sizeof(uint64_t));
	*attr = (uint64_t) big_ratio;
	trg->attr = attr;
	vcn_container_insert(big_trg, trg);
}

static inline msh_trg_t* refinement_remove_bigger_trg(vcn_container_t * restrict big_trg)
{
	int8_t status;
	msh_trg_t *trg = vcn_container_do(big_trg, "delete_last", NULL, &status);
	if (0 != status) {
		printf("\nError: vcn_mesh_refine()\n");
		printf("       in vcn_container_do(delete_last)\n");
		exit(1);
	}
	if (NULL != trg) {
		free(trg->attr);
		trg->attr = NULL;
	}
	return trg;  
}

static inline void refinement_remove_big_trg
                             (vcn_container_t *big_trg, msh_trg_t *trg)
{
	if (NULL != trg->attr) {
		if (NULL != vcn_container_delete(big_trg, trg)) {
			free(trg->attr);
			trg->attr = NULL;
		}
	}
}

static inline void refinement_check_trg
                          (msh_trg_t *const restrict trg,
			   vcn_mesh_t *const restrict mesh,
			   vcn_container_t *const restrict big_trg,
			   hash_trg_t *const restrict poor_quality_trg,
			   double cr2se_ratio,
			   double (*density)(const void *const data,
					     const double const x[2]),
			   const void *const restrict density_data)
{
	double big_ratio;
	if (mtrg_is_too_big(mesh, trg, density, density_data, &big_ratio)) {
		refinement_insert_big_trg(big_trg, trg, big_ratio);
	} else {
		if (cr2se_ratio > 1e3) 
			return;
		double trg_cr2se_ratio = vcn_utils2D_get_cr2se_ratio(trg->v1->x,
								  trg->v2->x,
								  trg->v3->x);
		if (trg_cr2se_ratio > cr2se_ratio)
			hash_trg_insert(poor_quality_trg, trg, trg_cr2se_ratio);
	}
}

static void refinement_initialize_big_and_poor_quality_trg
                          (vcn_mesh_t *const restrict mesh,
			   vcn_container_t *const restrict big_trg,
			   hash_trg_t *const restrict poor_quality_trg,
			   double cr2se_ratio,
			   double (*density)(const void *const data,
					     const double const x[2]),
			   const void *const restrict density_data)
{
	vcn_iterator_t *const restrict iter = vcn_iterator_create();
	vcn_iterator_set_container(iter, mesh->ht_trg);
	while (vcn_iterator_has_more(iter)) {
		msh_trg_t *restrict trg =
			(msh_trg_t*) vcn_iterator_get_next(iter);
		refinement_check_trg(trg, mesh, big_trg, poor_quality_trg,
				     cr2se_ratio, density, density_data);
	}
	vcn_iterator_destroy(iter);
}

static inline msh_trg_t* refinement_get_trg_containing_circumcenter
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

static inline vcn_container_t* refinement_get_sgm_encroached_by_vertex
                          (const vcn_mesh_t *const restrict mesh,
			   const msh_trg_t *const restrict trg_containing_vtx,
			   const msh_vtx_t *const restrict vtx)
{
	/* Get encroached triangles */
	vcn_container_t *const restrict encroached_trg =
		refinement_get_encroached_triangles(trg_containing_vtx, vtx);

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
	vcn_container_t *restrict encroached_sgm = vcn_container_create(VCN_CONTAINER_QUEUE);
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

static inline bool refinement_split_is_permitted
                          (const msh_edge_t *const restrict sgm,
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
		refinement_get_subsgm_cluster(sgm, sgm->v1, &smallest_angle);

	double smallest_angle_aux;
	vcn_container_t* cluster2 =
		refinement_get_subsgm_cluster(sgm, sgm->v2, &smallest_angle_aux);

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

static inline vcn_container_t* refinement_get_subsgm_cluster 
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

static void remove_triangles_propagate
                   (vcn_mesh_t *const restrict mesh, 
		    msh_trg_t* trg)
{
	mesh_substract_triangle(mesh, trg);
	trg->attr = (attr_t*)0x1;
	if (NULL != trg->t1) {
		msh_trg_t* nb_trg = trg->t1;
		mtrg_vanish_from_neighbour(trg, nb_trg);
		if (!medge_is_subsgm(trg->s1) && nb_trg->attr != (attr_t*)0x1)
			remove_triangles_propagate(mesh, nb_trg);
	}
	if (NULL != trg->t2) {
		msh_trg_t* nb_trg = trg->t2;
		mtrg_vanish_from_neighbour(trg, nb_trg);
		if (!medge_is_subsgm(trg->s2) && nb_trg->attr != (attr_t*)0x1)
			remove_triangles_propagate(mesh, nb_trg);
	}
	if (trg->t3 != NULL) {
		msh_trg_t* nb_trg = trg->t3;
		mtrg_vanish_from_neighbour(trg, nb_trg);
		if (!medge_is_subsgm(trg->s3) && nb_trg->attr != (attr_t*)0x1)
			remove_triangles_propagate(mesh, nb_trg);
	}
	free(trg);
}

static inline int compare_trg_attr_uint64_t(const void *const trgA,
					  const void *const trgB)
{
	register uint64_t a = *(uint64_t*)(((msh_trg_t*)trgA)->attr);
	register uint64_t b = *(uint64_t*)(((msh_trg_t*)trgB)->attr);
	if (a > b)
		return 1;
	else if (a < b)
		return -1;  
	return 0;
}

static inline int compare_deterministic_sgm
                 (const void *const restrict sgmA,
		  const void *const restrict sgmB)
{
	msh_edge_t *const restrict sA = (msh_edge_t*) sgmA;
	msh_edge_t *const restrict sB = (msh_edge_t*) sgmB;
	/* Compare x coordinate of vertex 1 */
	if (sA->v1->x[0] -  sB->v1->x[0] > VCN_GEOMETRIC_TOL)
		return 1;
	else if (sB->v1->x[0] - sA->v1->x[0] > VCN_GEOMETRIC_TOL)
		return -1;

	/* Compare y coordinate of vertex 1 */
	if (sA->v1->x[1] -  sB->v1->x[1] > VCN_GEOMETRIC_TOL)
		return 1;
	else if (sB->v1->x[1] - sA->v1->x[1] > VCN_GEOMETRIC_TOL)
		return -1;

	/* Compare x coordinate of vertex 2 */
	if (sA->v2->x[0] - sB->v2->x[0] > VCN_GEOMETRIC_TOL)
		return 1;
	else if (sB->v2->x[0] - sA->v2->x[0] > VCN_GEOMETRIC_TOL)
		return -1;

	/* Compare y coordinate of vertex 2 */
	if (sA->v2->x[1] - sB->v2->x[1] > VCN_GEOMETRIC_TOL)
		return 1;
	else if (sB->v2->x[1] - sA->v2->x[1] > VCN_GEOMETRIC_TOL)
		return -1;
  
	return 0;
}

static inline int compare_sgmA_isSmallerThan_sgmB
                                 (const void *const restrict sgmA, 
				  const void *const restrict sgmB)
{
	const uint64_t lA = 
		(uint64_t)(1e8 * medge_get_computed_length((msh_edge_t*)sgmA));
	const uint64_t lB = 
		(uint64_t)(1e8 * medge_get_computed_length((msh_edge_t*)sgmB));
	if (lB > lA)
		return 1;
	else if (lA > lB)
		return -1;
	return 0;
}

static inline int compare_sgmA_isTheSameThan_sgmB
                                 (const void *const restrict sgmA, 
				  const void *const restrict sgmB)
{
	const msh_edge_t *const restrict A = (msh_edge_t*)sgmA;
	const msh_edge_t *const restrict B = (msh_edge_t*)sgmB;
	if (A->v1 == B->v1 && A->v2 == B->v2)
		return 0;
	if (sgmA < sgmB)
		return 1;
	if (sgmA > sgmB) 
		return -1;
	return 0;
}

static inline int compare_trgA_isBetterThan_trgB
                                (const void *const restrict trgA,
				 const void *const restrict trgB)
{
	const uint64_t QA = 
		(uint64_t)(1e6 * mtrg_get_computed_quality((msh_trg_t*)trgA));
	const uint64_t QB = 
		(uint64_t)(1e6 * mtrg_get_computed_quality((msh_trg_t*)trgB));
	if (QA > QB)
		return 1;
	else if (QB > QA)
		return -1;  
	return 0;
}

static inline int compare_trgA_isSmallerThan_trgB
                                (const void *const restrict trgA, 
				 const void *const restrict trgB)
{
	const uint64_t SA = 
		(uint64_t)(1e8 * mtrg_get_computed_size((msh_trg_t*)trgA));
	const uint64_t SB = 
		(uint64_t)(1e8 * mtrg_get_computed_size((msh_trg_t*)trgB));
	if (SA < SB)
		return 1;
	else if (SB < SA)
		return -1;  
	return 0;
}

static inline int compare_area1_isGreaterThan_area2
                                (const void *const restrict a1,
				 const void *const restrict a2)
{
	double area1_d = ((double*)((void**)a1)[0])[0];
	double area2_d = ((double*)((void**)a2)[0])[0];
	uint64_t area1 = (uint64_t)(1e8 * area1_d);
	uint64_t area2 = (uint64_t)(1e8 * area2_d);
	if (area2 < area1) 
		return 1;
	else if (area2 > area1)
		return -1;
	else
		return 0;
}

static inline int compare_id(const void* const ptrA,
			     const void* const ptrB)
{
	uint32_t* A = (uint32_t*) ptrA;
	uint32_t* B = (uint32_t*) ptrB;
	if (A[0] == B[0])
		return 0;
	if(A[0] > B[0])
		return 1;
	return -1;
}

static inline int compare_vtx(const void* const vtxA,
			      const void* const vtxB)
{
	msh_vtx_t* vA = (msh_vtx_t*) vtxA;
	msh_vtx_t* vB = (msh_vtx_t*) vtxB;
	if (vcn_utils2D_get_dist2(vA->x, vB->x) < VCN_GEOMETRIC_TOL)
		return 0;
	return 1;
}

static inline int compare_sgm_by_dist_and_by_vtx_pointer
                                            (const void* const sgmA,
					     const void* const sgmB)
{
	int comp1 = compare_sgmA_isSmallerThan_sgmB(sgmA, sgmB);
	if (0 == comp1) {
		return compare_sgmA_isTheSameThan_sgmB(sgmA, sgmB);
	}
	return -comp1;
}

static inline uint32_t hash_key_vtx(const void *const  vertex)
{
	const msh_vtx_t *const vtx = vertex;
	return (uint32_t)((int)(vtx->x[0] * 73856093) ^
			  (int)(vtx->x[1] * 19349663));
}

static inline uint32_t hash_key_trg(const void *const restrict triangle)
{
	const msh_trg_t *const restrict trg = triangle;
	return (uint32_t)
		((int)((trg->v1->x[0]/((trg->v1->x[1]!=0.0)?
				       trg->v1->x[1]:13.13))*73856093) ^ 
		 (int)((trg->v2->x[1]/((trg->v2->x[0]!=0.0)?
				       trg->v2->x[0]:17.17))*19349663) ^
		 (int)((trg->v3->x[0]/((trg->v3->x[1]!=0.0)?
				       trg->v3->x[1]:23.23))*83492791));
}

static void remove_concavities_triangles(vcn_mesh_t* mesh)
{
	/* (BIG OPPORTUNITY TO MAKE IT FAST) */
	msh_trg_t* trg = NULL;
	vcn_iterator_t* iter = vcn_iterator_create();
	vcn_iterator_set_container(iter, mesh->ht_trg);
	while (vcn_iterator_has_more(iter)) {
		trg = (msh_trg_t*)vcn_iterator_get_next(iter);
		if ((trg->t1 == NULL && !medge_is_subsgm(trg->s1)) || 
		    (trg->t2 == NULL && !medge_is_subsgm(trg->s2)) ||
		    (trg->t3 == NULL && !medge_is_subsgm(trg->s3))) {
			remove_triangles_propagate(mesh, trg);
			vcn_iterator_restart(iter);
		}
	}
	vcn_iterator_destroy(iter);
}

static void remove_holes_triangles(vcn_mesh_t* mesh, double* holes, uint32_t N_holes)
{
	/* (BIG OPPORTUNITY TO MAKE IT FAST) */
	for (uint32_t i=0; i < N_holes; i++) {
		msh_trg_t* trg = NULL;
		vcn_iterator_t* iter = vcn_iterator_create();
		vcn_iterator_set_container(iter, mesh->ht_trg);
		while (vcn_iterator_has_more(iter)) {
			/* REFACTOR next line */
			trg = (msh_trg_t*) vcn_iterator_get_next(iter);
			if (vcn_utils2D_pnt_lies_in_trg(trg->v1->x,
							trg->v2->x,
							trg->v3->x,
							&(holes[i*2])))
			{
				remove_triangles_propagate(mesh, trg);
				break;
			}
		}
		vcn_iterator_destroy(iter);
	}
}

void vcn_mesh_refine(vcn_mesh_t *const restrict mesh,
		     uint32_t max_vtx,
		     uint32_t max_trg,
		     double min_angle,
		     double (*density)(const void *const data,
				       const double const x[2]),
		     const void *const density_data)
{
	/* Allocate data structures to allocate encroached elements */
	vcn_container_t *restrict encroached_sgm =
		vcn_container_create(VCN_CONTAINER_SORTED);
	/* REFACTOR
	   set_key_generator(compare_deterministic_sgm);
	   set_comparer(...)
	*/
	vcn_container_t *restrict big_trg =
		vcn_container_create(VCN_CONTAINER_SORTED);
	/* REFACTOR
	   set_key_generator(compare_trg_attr_uint64_t);
	   set_comparer(...)
	*/
	hash_trg_t *restrict poor_quality_trg = hash_trg_create();

	/* Initialize FIFO with encroached segments */
	refinement_initialize_encroached_sgm(mesh, encroached_sgm);

	/* Calculate max circumradius to shortest edge ratio allowed */
	if (min_angle > VCN_ANGLE_MAX) {
		if (max_vtx == 0) {
			printf("WARNING in vcn_mesh_refine(): Setting max_vtx = 1000000 as\
              guarantee to finish.\n");
			/* Set a maximum because there are not guarantees to finish */
			max_vtx = 1000000;
		}
	}
	double cr2se_ratio = 1e30;
	if (min_angle > 0.0)
		cr2se_ratio = 1.0 / (2.0 * sin(min_angle));

	/* Split encroached segments */
	refinement_split_encroached_segments(mesh, encroached_sgm,
					     big_trg, poor_quality_trg,
					     1e30, NULL, NULL);

	/* Initialize hash table with the poor quality triangles */
	refinement_initialize_big_and_poor_quality_trg(mesh, big_trg,
						       poor_quality_trg,
						       cr2se_ratio,
						       density, density_data);

	/* Remove poor quality triangles */
	uint32_t iter = 0;
	while ((vcn_container_is_not_empty(big_trg) ||
		hash_trg_length(poor_quality_trg) > 0) &&
	       refinement_check_max_trg(mesh, max_trg) &&
	       refinement_check_max_vtx(mesh, max_vtx)) {
		/* Re-allocate hash tables if they are too small */
		if (iter % 5000 == 0) {
			if (vcn_bins2D_get_min_points_x_bin(mesh->ug_vtx) > 100)
				mesh_reallocate_htables(mesh);
		}
		iter ++;
		/* Get next poor-quality triangle */
		msh_trg_t *trg;
		if (vcn_container_is_not_empty(big_trg))
			trg = refinement_remove_bigger_trg(big_trg);
		else
			trg = hash_trg_remove_first(poor_quality_trg);

		/* Get circumcenter */
		msh_vtx_t *restrict cc = (msh_vtx_t*) vcn_point2D_create();
		vcn_utils2D_get_circumcenter(trg->v1->x, trg->v2->x, trg->v3->x, cc->x);
    
		msh_trg_t *restrict trg_containing_cc =
			refinement_get_trg_containing_circumcenter(mesh, trg, cc);

		vcn_container_t *restrict  sgm_encroached_by_cc = 
			refinement_get_sgm_encroached_by_vertex(mesh, trg_containing_cc, cc);
		if (vcn_container_is_empty(sgm_encroached_by_cc)) {
			/* Circumcenter accepted */
			vcn_container_t *restrict l_new_trg = 
				vcn_container_create(VCN_CONTAINER_QUEUE);

			refinement_insert_circumcenter(mesh, trg_containing_cc, cc, big_trg,
						       poor_quality_trg, l_new_trg);

			refinement_verify_new_encroachments(mesh, cc, l_new_trg, NULL, big_trg,
							    poor_quality_trg, cr2se_ratio,
							    density, density_data);
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
				if (mtrg_is_too_big(mesh, trg, density, density_data, NULL))
					vcn_container_insert(encroached_sgm, sgm);
				else if (refinement_split_is_permitted(sgm, d))
					vcn_container_insert(encroached_sgm, sgm);
			}
			/* Process encroached segments */
			if (vcn_container_is_not_empty(encroached_sgm)) {
				/* Put it back for another try */
				refinement_check_trg(trg, mesh, big_trg, poor_quality_trg,
						     cr2se_ratio, density, density_data);
	
				/* Split encroached segments */
				refinement_split_encroached_segments(mesh, encroached_sgm,
								     big_trg, poor_quality_trg,
								     cr2se_ratio,
								     density,
								     density_data);
			}
		}
		vcn_container_destroy(sgm_encroached_by_cc);
	}
	/* Free data structures */
	vcn_container_destroy(encroached_sgm);
	vcn_container_destroy(big_trg);
	hash_trg_destroy(poor_quality_trg);
}

int vcn_mesh_insert_vertex(vcn_mesh_t *const restrict mesh, 
			   const double *const restrict vertex)
{
	msh_vtx_t* new_vtx = (msh_vtx_t*) vcn_point2D_create();

	/* Move and scale new vertex */
	new_vtx->x[0] = mesh->scale * (vertex[0] - mesh->xdisp);
	new_vtx->x[1] = mesh->scale * (vertex[1] - mesh->ydisp);
  
	msh_trg_t* trg = mesh_locate_vtx(mesh, new_vtx);

	if (NULL == trg) {
		vcn_point2D_destroy(new_vtx);
		return 1;
	}

	refinement_insert_circumcenter(mesh, trg, new_vtx, NULL, NULL, NULL);
	return 0;
}

vcn_mesh_t* vcn_mesh_clone(const vcn_mesh_t* const mesh)
{
	vcn_mesh_t* clone = malloc(sizeof(*clone));

	/* Copy scale and shift values to handle floating point error */
	clone->scale = mesh->scale;
	clone->xdisp = mesh->xdisp;
	clone->ydisp = mesh->ydisp;

	/* Clone arrays to relate the mesh with the input */
	clone->N_input_vtx = mesh->N_input_vtx;
	clone->input_vtx = 
		calloc(clone->N_input_vtx, sizeof(*(clone->input_vtx)));
	clone->N_input_sgm = mesh->N_input_sgm;
	clone->input_sgm =
		calloc(clone->N_input_sgm, sizeof(*(clone->input_sgm)));

	/* Clone grid of vertices */
	uint32_t N_vertices = vcn_bins2D_get_length(mesh->ug_vtx);
	msh_vtx_t** vertices = malloc(N_vertices * sizeof(*vertices));
	double bins_size = vcn_bins2D_get_size_of_bins(mesh->ug_vtx);
	clone->ug_vtx = vcn_bins2D_create(bins_size);
	uint32_t i = 0;
	vcn_bins2D_iter_t* giter = vcn_bins2D_iter_create();
	vcn_bins2D_iter_set_bins(giter, mesh->ug_vtx);
	while (vcn_bins2D_iter_has_more(giter)) {
		msh_vtx_t* vtx = (msh_vtx_t*) vcn_bins2D_iter_get_next(giter);
		/* Create the vertex clone */
		msh_vtx_t* vtx_clone = vcn_point2D_create();
		memcpy(vtx_clone->x, vtx->x, 2 * sizeof(*(vtx->x)));
		vtx_clone->attr = vtx->attr;

		/* Set an index to the original vertex */
		void** attr = malloc(2 * sizeof(*attr));
		uint32_t *id = malloc(sizeof(*id));
		id[0] = i++;
		attr[0] = id;
		attr[1] = vtx->attr;
		vtx->attr = attr;
		/* Set cloned and original vertices to the built-in hash table */
		vertices[id[0]] = vtx_clone;
		/* Insert the clone vertex into the clone grid */
		vcn_bins2D_insert(clone->ug_vtx, vtx_clone);
	}
	/* Create a built-in hash table to relate original and cloned segments
	 * and triangles */
	uint32_t N_triangles = vcn_container_get_length(mesh->ht_trg);
	msh_trg_t** triangles = malloc(N_triangles * sizeof(*triangles));
	i = 0;
	vcn_iterator_t* trg_iter = vcn_iterator_create();
	vcn_iterator_set_container(trg_iter, mesh->ht_trg);
	while (vcn_iterator_has_more(trg_iter)) {
		msh_trg_t* trg = (msh_trg_t*) vcn_iterator_get_next(trg_iter);
		/* Clone the triangle */
		msh_trg_t* trg_clone = calloc(1, sizeof(*trg_clone));
		trg_clone->attr = trg->attr;
		/* Set an index to the original */
		void** attr = malloc(2 * sizeof(*attr));
		uint32_t *id = malloc(sizeof(*id));
		id[0] = i++;
		attr[0] = id;
		attr[1] = trg->attr;
		trg->attr = attr;
		/* Set cloned and original to the built-in hash table */
		triangles[id[0]] = trg_clone;
	}

	uint32_t N_segments = vcn_container_get_length(mesh->ht_edge);
	msh_edge_t** segments = malloc(N_segments * sizeof(*segments));
	i = 0;
	vcn_iterator_t* sgm_iter = vcn_iterator_create();
	vcn_iterator_set_container(sgm_iter, mesh->ht_edge);
	while (vcn_iterator_has_more(sgm_iter)) {
		msh_edge_t* sgm = (msh_edge_t*) vcn_iterator_get_next(sgm_iter);
		/* Clone the segment */
		msh_edge_t* sgm_clone = calloc(1, sizeof(*sgm_clone));

		/* Set an index to the original */
		void** attr = malloc(2 * sizeof(*attr));
		uint32_t *id = malloc(sizeof(*id));
		id[0] = i++;
		attr[0] = id;
		attr[1] = sgm->attr;
		sgm->attr = attr;
		/* Set cloned and original to the built-in hash table */
		segments[id[0]] = sgm_clone;
	}

	/* Clone data structures and hash tables used to handle mesh */
	clone->ht_trg = vcn_container_create(VCN_CONTAINER_HASH);
	vcn_container_set_key_generator(clone->ht_trg, hash_key_trg);
	vcn_iterator_restart(trg_iter);
	while (vcn_iterator_has_more(trg_iter)) {
		const msh_trg_t *trg = vcn_iterator_get_next(trg_iter);
		msh_trg_t *trg_clone = 
			triangles[((uint32_t*)((void**)trg->attr)[0])[0]];
		trg_clone->v1 = vertices[((uint32_t*)((void**)trg->v1->attr)[0])[0]];
		trg_clone->v2 = vertices[((uint32_t*)((void**)trg->v2->attr)[0])[0]];
		trg_clone->v3 = vertices[((uint32_t*)((void**)trg->v3->attr)[0])[0]];

		trg_clone->s1 = segments[((uint32_t*)((void**)trg->s1->attr)[0])[0]];
		trg_clone->s2 = segments[((uint32_t*)((void**)trg->s2->attr)[0])[0]];
		trg_clone->s3 = segments[((uint32_t*)((void**)trg->s3->attr)[0])[0]];

		if (NULL != trg->t1)
			trg_clone->t1 = triangles[((uint32_t*)((void**)trg->t1->attr)[0])[0]];
		if (NULL != trg->t2)
			trg_clone->t2 = triangles[((uint32_t*)((void**)trg->t2->attr)[0])[0]];
		if (NULL != trg->t3)
			trg_clone->t3 = triangles[((uint32_t*)((void**)trg->t3->attr)[0])[0]];
		vcn_container_insert(clone->ht_trg, trg_clone);
	}
	clone->ht_edge = vcn_container_create(VCN_CONTAINER_HASH);
	vcn_container_set_key_generator(clone->ht_edge, hash_key_edge);
	vcn_container_set_comparer(clone->ht_edge, are_equal_edge);
	vcn_iterator_restart(sgm_iter);
	while (vcn_iterator_has_more(sgm_iter)) {
		const msh_edge_t *sgm = vcn_iterator_get_next(sgm_iter);
		msh_edge_t *sgm_clone = segments[((uint32_t*)((void**)sgm->attr)[0])[0]];
		sgm_clone->v1 = vertices[((uint32_t*)((void**)sgm->v1->attr)[0])[0]];
		sgm_clone->v2 = vertices[((uint32_t*)((void**)sgm->v2->attr)[0])[0]];
		if (NULL != sgm->t1)
			sgm_clone->t1 = triangles[((uint32_t*)((void**)sgm->t1->attr)[0])[0]];
		if (NULL != sgm->t2)
			sgm_clone->t2 = triangles[((uint32_t*)((void**)sgm->t2->attr)[0])[0]];
		attr_t* sgm_attr = (attr_t*)((void**)sgm->attr)[1];
		if (NULL != sgm_attr) {
			if (1 == sgm_attr->id) {
				/* The segment is forming the input */
				input_sgm_attr_t* input_sgm_attr = sgm_attr->data;
				msh_edge_t* prev_clone = NULL;
				msh_edge_t* prev = input_sgm_attr->prev;
				if (NULL != prev)
					prev_clone = segments[((uint32_t*)((void**)prev->attr)[0])[0]];
				msh_edge_t* next_clone = NULL;
				msh_edge_t* next = input_sgm_attr->next;
				if(NULL != next)
					next_clone = segments[((uint32_t*)((void**)next->attr)[0])[0]];
				medge_set_as_subsgm(sgm_clone, medge_subsgm_idx(sgm),
						    prev_clone, next_clone);
			}
		}
		vcn_container_insert(clone->ht_edge, sgm_clone);
	}
	for (uint32_t i = 0; i < clone->N_input_vtx; i++)
		clone->input_vtx[i] = 
			vertices[((uint32_t*)((void**)mesh->input_vtx[i]->attr)[0])[0]];

	for (uint32_t i = 0; i < clone->N_input_sgm; i++)
		clone->input_sgm[i] = 
			segments[((uint32_t*)((void**)mesh->input_sgm[i]->attr)[0])[0]];  

	/* Free memory */
	free(vertices);
	free(segments);
	free(triangles);

	vcn_iterator_restart(sgm_iter);
	while (vcn_iterator_has_more(sgm_iter)) {
		msh_edge_t* sgm = (msh_edge_t*) vcn_iterator_get_next(sgm_iter);
		void** attr = sgm->attr;
		sgm->attr = attr[1];
		free(attr[0]);
		free(attr);
	}
	vcn_iterator_destroy(sgm_iter);

	vcn_iterator_restart(trg_iter);
	while (vcn_iterator_has_more(trg_iter)) {
		msh_trg_t* trg = (msh_trg_t*) vcn_iterator_get_next(trg_iter);
		void** attr = trg->attr;
		trg->attr = attr[1];
		free(attr[0]);
		free(attr);
	}
	vcn_iterator_destroy(trg_iter);

	vcn_bins2D_iter_restart(giter);
	while (vcn_bins2D_iter_has_more(giter)) {
		msh_vtx_t* vtx = (msh_vtx_t*) vcn_bins2D_iter_get_next(giter);
		void** attr = vtx->attr;
		vtx->attr = attr[1];
		free(attr[0]);
		free(attr);
	}
	vcn_bins2D_iter_destroy(giter);

	/* Return clone */
	return clone;
}

void vcn_mesh_generate_from_model(vcn_mesh_t *mesh,
				  const vcn_model_t *const restrict model)
{
	;/* REFACTOR */
}

vcn_mesh_t* vcn_mesh_create_from_model(const vcn_model_t *const restrict model,
				       uint32_t max_vtx,
				       uint32_t max_trg,
				       double min_angle,
				       double (*density)(const void *const data,
							 const double const x[2]),
				       const void *const restrict density_data)
{
	/* Create Constrained Delaunay triangulation */
	vcn_mesh_t* mesh = vcn_mesh_get_constrained_delaunay(model->N,
							     model->vertex,
							     model->M,
							     model->edge);

	/* Allocate scalled holes */
	double* holes = NULL; /* Dummy initialization */
	if (0 < model->H) {
		holes = calloc(2 * model->H, sizeof(*holes));
		memcpy(holes, model->holes, 2 * model->H * sizeof(*holes));
	}

	/* Removing holes and concavities */
	for (uint32_t i = 0; i < model->H; i++) {
		holes[i * 2] = mesh->scale * (model->holes[i * 2] - mesh->xdisp);
		holes[i*2+1] = mesh->scale * (model->holes[i*2+1] - mesh->ydisp);
	}

	remove_concavities_triangles(mesh);
	if (0 < model->H) {
		remove_holes_triangles(mesh, holes, model->H);
		free(holes);
	}

	/* Refine mesh */
	if (density != VCN_DENSITY_CDT) {
		if (vcn_bins2D_get_length(mesh->ug_vtx) < max_vtx || max_vtx == 0)
			if (vcn_container_get_length(mesh->ht_trg) < max_trg || max_trg == 0)
				vcn_mesh_refine(mesh, max_vtx, max_trg, min_angle,
						density, density_data);
	}
	/* Return mesh */
	return mesh;
}
