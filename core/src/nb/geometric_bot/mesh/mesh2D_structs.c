#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <alloca.h>

#include "nb/math_bot.h"
#include "nb/container_bot/container.h"
#include "nb/container_bot/iterator.h"
#include "nb/geometric_bot/utils2D.h"
#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/geometric_bot/knn/bins2D.h"
#include "nb/geometric_bot/knn/bins2D_iterator.h"

#include "mesh2D_structs.h"

static bool mesh_remove_edge
                       (nb_container_t *const ht_edge,
			const msh_vtx_t *const v1,
			const msh_vtx_t *const v2);


bool medge_is_boundary(const msh_edge_t *const edge)
{
	return (NULL == edge->t1 || NULL == edge->t2);
}

inline bool medge_is_subsgm(const msh_edge_t *const restrict sgm)
{
	if (NULL == sgm->attr)
		return false;
	return ((attr_t*)sgm->attr)->id == 1;
}

void medge_set_as_subsgm(msh_edge_t *const restrict sgm,
			 uint32_t idx,
			 const msh_edge_t *const restrict prev,
			 const msh_edge_t *const restrict next)
{
	attr_t *const restrict sgm_attr = (attr_t*)calloc(1, sizeof(attr_t));
	sgm_attr->id = 1;
	input_sgm_attr_t *const restrict attr = 
		(input_sgm_attr_t*)calloc(1, sizeof(input_sgm_attr_t));
	attr->idx = idx;
	attr->prev = (msh_edge_t*)prev;
	attr->next = (msh_edge_t*)next;
	sgm_attr->data = attr;
	attr->attr = sgm->attr;
	sgm->attr = sgm_attr;
}

inline void medge_update_subsgm_next(msh_edge_t *const restrict sgm,
				     const msh_edge_t *const restrict next)
{
	input_sgm_attr_t *const restrict attr = 
		(input_sgm_attr_t*)((attr_t*)sgm->attr)->data; 
	attr->next = (msh_edge_t*)next;
}

inline void medge_update_subsgm_prev(msh_edge_t *const sgm, 
				     const msh_edge_t *const prev)
{
	input_sgm_attr_t* attr = 
		(input_sgm_attr_t*)((attr_t*)sgm->attr)->data;
	attr->prev = (msh_edge_t*)prev; 
	/* Casting from const to non-const pointer
	 * in order to link it to the double linked
	 * list, but with the compromise to do not 
	 * modify its value */
}

inline uint32_t medge_subsgm_idx(const msh_edge_t *const restrict sgm)
{
	input_sgm_attr_t* attr = 
		(input_sgm_attr_t*)((attr_t*)sgm->attr)->data;
	return attr->idx;
}

inline msh_edge_t* medge_subsgm_prev(const msh_edge_t *const restrict sgm)
{
	input_sgm_attr_t* attr = 
		(input_sgm_attr_t*)((attr_t*)sgm->attr)->data;
	return attr->prev;
}

inline msh_edge_t* medge_subsgm_next(const msh_edge_t *const restrict sgm)
{
	input_sgm_attr_t* attr = 
		(input_sgm_attr_t*)((attr_t*)sgm->attr)->data;
	return attr->next;
}

void medge_subsgm_get_extreme_vtx(const msh_edge_t *const sgm,
				  msh_vtx_t *vtx[2])
{
	msh_edge_t *first = sgm;
	while (NULL != medge_subsgm_prev(first))
		first = medge_subsgm_prev(first);

	msh_edge_t *last = sgm;
	while (NULL != medge_subsgm_next(last))
		last = medge_subsgm_next(last);

	if (first == last) {
		vtx[0] = sgm->v1;
		vtx[1] = sgm->v2;
	} else {
		msh_edge_t *fnext = medge_subsgm_next(first);
		vtx[0] =  medge_subsgm_get_opposite_vtx(first, fnext);

		msh_edge_t *lprev = medge_subsgm_prev(last);
		vtx[1] =  medge_subsgm_get_opposite_vtx(last, lprev);
	}
}

msh_vtx_t* medge_subsgm_get_opposite_vtx(const msh_edge_t *const sgm,
					 const msh_edge_t *const adj)
{
	msh_vtx_t *vtx;
	if (sgm->v1 == adj->v1 || sgm->v1 == adj->v2)
		vtx = sgm->v2;
	else if (sgm->v2 == adj->v1 || sgm->v2 == adj->v2)
		vtx = sgm->v1;
	else
		vtx = NULL;
	return vtx;
}

inline void medge_subsgm_set_attribute(msh_edge_t* sgm, void* attr)
{
	input_sgm_attr_t* in_attr = 
		(input_sgm_attr_t*)((attr_t*)sgm->attr)->data;
	in_attr->attr = attr;
}

inline void* medge_subsgm_get_attribute(const msh_edge_t* const sgm)
{
	input_sgm_attr_t* in_attr = 
		(input_sgm_attr_t*)((attr_t*)sgm->attr)->data;
	return in_attr->attr;
}

void medge_destroy_subsgm_attribute(msh_edge_t *const restrict sgm)
{
	void* attr = ((input_sgm_attr_t*)((attr_t*)sgm->attr)->data)->attr;
	free(((attr_t*)sgm->attr)->data);
	free(sgm->attr);
	sgm->attr = attr;
}

void medge_set_length(msh_edge_t *const sgm)
{
	attr_t* attr = (attr_t*)calloc(1, sizeof(attr_t));
	attr->id = 2;
	double* length = malloc(sizeof(*length));
	*length = vcn_utils2D_get_dist(sgm->v1->x, sgm->v2->x);
	attr->data = length;

	if(medge_is_subsgm(sgm))
		medge_subsgm_set_attribute(sgm, attr);
	else
		sgm->attr = attr;
}

void medge_update_length(msh_edge_t *const sgm)
{
	attr_t* attr;
	if (medge_is_subsgm(sgm))
		attr = (attr_t*) medge_subsgm_get_attribute(sgm);
	else
		attr = (attr_t*) sgm->attr;
	if (NULL == attr)
		/* TEMPORAL: This should never happen
		 * The routine mesh_insert_vtx does not allocate length 
		 */
		return;
	if (attr->id != 2) 
		return;
	double* length = attr->data;
	length[0] = vcn_utils2D_get_dist(sgm->v1->x, sgm->v2->x);
}

double medge_get_computed_length(const msh_edge_t *const sgm)
{
	attr_t* attr;
	if (medge_is_subsgm(sgm))
		attr = (attr_t*) medge_subsgm_get_attribute(sgm);
	else
		attr = (attr_t*) sgm->attr;
	if (NULL == attr)
		/* TEMPORAL: This should never happen
		 * The routine mesh_insert_vtx does not allocate length 
		 */
		return vcn_utils2D_get_dist(sgm->v1->x, sgm->v2->x);
	if (attr->id != 2)
		return vcn_utils2D_get_dist(sgm->v1->x, sgm->v2->x);
	double* length = attr->data;
	return length[0];
  
}

void medge_destroy_length_attribute(msh_edge_t *const sgm)
{
	attr_t* attr;
	if (medge_is_subsgm(sgm))
		attr = medge_subsgm_get_attribute(sgm);
	else
		attr = sgm->attr;

	if (NULL == attr)
		/* TEMPORAL: This should never happen */
		return;
	if (attr->id != 2) 
		return;

	if (medge_is_subsgm(sgm))
		medge_subsgm_set_attribute(sgm, NULL);
	else
		sgm->attr = NULL;

	free(attr->data);
	free(attr);
}

bool medge_has_a_triangle_to_the_left
                             (const msh_edge_t *const restrict sgm, 
			      const msh_vtx_t *const restrict v1)
{
	/* Important: The vertex is just to know the direction */
	if (NULL == sgm)
		return false;
	if (sgm->v1 == v1) {
		if (sgm->t1 != NULL)
			return true;
		else return false;
	} else /* sgm->v2 == v1 */{
		if (sgm->t2 != NULL)
			return true;
		else return false;
	}
}

void medge_connect_triangles(msh_edge_t *const restrict sgm)
{
	if (sgm->t1 != NULL && sgm->t2 != NULL) {
		/* Create triangles connection */
		msh_trg_t* t1 = sgm->t1;
		msh_trg_t* t2 = sgm->t2;
		if (sgm->v1 == t1->v1)
			t1->t1 = t2;
		else if (sgm->v1 == t1->v2)
			t1->t2 = t2;
		else if (sgm->v1 == t1->v3)
			t1->t3 = t2;

		if (sgm->v2 == t2->v1)
			t2->t1 = t1;
		else if (sgm->v2 == t2->v2)
			t2->t2 = t1;
		else if (sgm->v2 == t2->v3)
			t2->t3 = t1;
	}
}

void medge_connect_triangle
                         (msh_edge_t *const restrict sgm, 
			  msh_trg_t *const trg,
			  const msh_vtx_t *const restrict v1, 
			  const msh_vtx_t *const restrict v2)
{
	if (sgm->v1 == v1) {
		sgm->t1 = trg;
		if (sgm->v1 == trg->v1)
			trg->s1 = sgm;
		else if (sgm->v1 == trg->v2)
			trg->s2 = sgm;
		else if (sgm->v1 == trg->v3)
			trg->s3 = sgm;
	} else if (sgm->v1 == v2) {
		sgm->t2 = trg;
		if (sgm->v1 == trg->v2)
			trg->s1 = sgm;
		else if (sgm->v1 == trg->v3)
			trg->s2 = sgm;
		else if (sgm->v1 == trg->v1)
			trg->s3 = sgm;
	}
  
	if (sgm->t1 != NULL && sgm->t2 != NULL) {
		/* Create triangles connection */
		msh_trg_t* t1 = sgm->t1;
		msh_trg_t* t2 = sgm->t2;
		if (sgm->v1 == t1->v1)
			t1->t1 = t2;
		else if (sgm->v1 == t1->v2)
			t1->t2 = t2;
		else if (sgm->v1 == t1->v3)
			t1->t3 = t2;

		if (sgm->v2 == t2->v1)
			t2->t1 = t1;
		else if (sgm->v2 == t2->v2)
			t2->t2 = t1;
		else if (sgm->v2 == t2->v3)
			t2->t3 = t1;
	}
}

msh_trg_t* medge_get_triangle_to_the_left
                           (const msh_edge_t *const restrict sgm,
			    const msh_vtx_t *const restrict v1)
{
	/* Important: The vertex is just to know the direction */
	if(NULL == sgm)
		return NULL;
	if(sgm->v1 == v1)
		return (msh_trg_t*)sgm->t1;
	if(sgm->v2 == v1)
		return (msh_trg_t*)sgm->t2;
	return NULL;
}

inline msh_trg_t* medge_get_opposite_triangle
                                  (const msh_edge_t *const sgm,
				   const msh_trg_t *const trg)
{
	if (trg == sgm->t1)
		return sgm->t2;
	if (trg == sgm->t2)
		return sgm->t1;
	return NULL;
}

inline msh_trg_t* mtrg_create(void)
{
	return (msh_trg_t*) calloc(1, sizeof(msh_trg_t));
}

void mtrg_set_quality_and_size(msh_trg_t *const restrict trg)
/* This function assumes that the segments has the length stored as an attribute */
{
  /* Set attribute */
  attr_t* attr = (attr_t*) calloc(1, sizeof(attr_t));
  attr->id = 1;
  double* data = (double*) malloc(2 * sizeof(double));

  /* Compute the double of the area */
  const double area2 = vcn_utils2D_get_2x_trg_area(trg->v1->x, trg->v2->x, trg->v3->x);

  /* Calculate quality */
  const double circumradius = 
    (medge_get_computed_length(trg->s1) *
     medge_get_computed_length(trg->s2) *
     medge_get_computed_length(trg->s3))/(2.0 * area2);
  const msh_edge_t* min_sgm = mtrg_get_shortest_edge(trg);
  const double min_length = medge_get_computed_length(min_sgm);
  data[0] = NB_MATH_INV_SQRT3 * (min_length/circumradius);

  /* Calculate area */
  data[1] = 0.5 * area2;

  attr->data = data;
  trg->attr = attr;
}

double mtrg_get_computed_quality(const msh_trg_t *const restrict trg)
{
  if(trg->attr == NULL) 
    /* TEMPORAL: This should never happen
     * The func mesh_insert_vertex doesn't allocate quality
     */
    return vcn_utils2D_get_trg_quality(trg->v1->x,
				       trg->v2->x,
				       trg->v3->x);
  attr_t* attr = (attr_t*) trg->attr;
  if(attr->id != 1) 
    /* TEMPORAL: This should never happen
     * The func mesh_insert_vertex doesn't allocate quality
     */
    return vcn_utils2D_get_trg_quality(trg->v1->x,
				       trg->v2->x,
				       trg->v3->x);
  return ((double*)attr->data)[0];
}

double mtrg_get_computed_size(const msh_trg_t *const restrict trg)
{
  if(trg->attr == NULL)
    /* TEMPORAL: This should never happen
     * The func mesh_insert_vertex doesn't allocate quality
     */
    return vcn_utils2D_get_2x_trg_area(trg->v1->x,
				       trg->v2->x,
				       trg->v3->x);
  attr_t* attr = (attr_t*) trg->attr;
  if(attr->id != 1)
    /* TEMPORAL: This should never happen
     * The func mesh_insert_vertex doesn't allocate quality
     */
    return vcn_utils2D_get_2x_trg_area(trg->v1->x,
				       trg->v2->x,
				       trg->v3->x);
  return ((double*)attr->data)[1];
}

inline void mtrg_destroy_quality_and_size_attributes
                  (msh_trg_t *const trg)
{
	if (NULL == trg->attr)
		return;
	attr_t* attr = (attr_t*) trg->attr;
	if (attr->id != 1)
		return;
	free(attr->data);
	free(attr);
	trg->attr = NULL;
}

inline bool mtrg_has_an_input_vertex(const msh_trg_t *const trg)
{
	if (trg->v1->attr == (void*)1)
		return true;
	if (trg->v2->attr == (void*)1)
		return true;
	if (trg->v3->attr == (void*)1)
		return true;
	return false;
}

inline msh_edge_t* mtrg_get_largest_edge(const msh_trg_t *const trg)
{
	/* This function assumes that the length is stored as an attribute */
	msh_edge_t* sgm = trg->s1;
	const double l1 = medge_get_computed_length(trg->s1);
	const double l2 = medge_get_computed_length(trg->s2);
	const double l3 = medge_get_computed_length(trg->s3);
	if (l2 > l3) {
		if (l2 > l1)
			sgm = trg->s2;
	} else {
		if (l3 > l1)
			sgm = trg->s3;
	}    
	return sgm;
}

inline msh_edge_t* mtrg_get_shortest_edge(const msh_trg_t *const trg)
{
	/* This function assumes that the length is stored as an attribute */
	msh_edge_t* sgm = trg->s1;
	const double l1 = medge_get_computed_length(trg->s1);
	const double l2 = medge_get_computed_length(trg->s2);
	const double l3 = medge_get_computed_length(trg->s3);
	if (l2 < l3) {
		if (l2 < l1)
			sgm = trg->s2;
	} else {
		if (l3 < l1)
			sgm = trg->s3;
	}    
	return sgm;
}

msh_vtx_t* mtrg_get_opposite_vertex_guided
                           (const msh_trg_t *const restrict trg, 
			    const msh_edge_t *const restrict sgm,
			    bool same_vertices_order)
{
	if (same_vertices_order) {
		if (trg->v1 == sgm->v1)
			return trg->v3;
		if (trg->v2 == sgm->v1)
			return trg->v1;
		if (trg->v3 == sgm->v1)
			return trg->v2;
	} else {
		if (trg->v1 == sgm->v2)
			return trg->v3;
		if (trg->v2 == sgm->v2)
			return trg->v1;
		if (trg->v3 == sgm->v2)
			return trg->v2;
	}
	return NULL;
}

msh_vtx_t* mtrg_get_opposite_vertex(const msh_trg_t *const restrict trg, 
				    const msh_edge_t *const restrict sgm)
{
	if (trg->v1 == sgm->v1 && trg->v2 == sgm->v2)
		return trg->v3;
	if (trg->v2 == sgm->v1 && trg->v3 == sgm->v2)
		return trg->v1;
	if (trg->v3 == sgm->v1 && trg->v1 == sgm->v2)
		return trg->v2;
	if (trg->v1 == sgm->v2 && trg->v2 == sgm->v1)
		return trg->v3;
	if (trg->v2 == sgm->v2 && trg->v3 == sgm->v1)
		return trg->v1;
 	if (trg->v3 == sgm->v2 && trg->v1 == sgm->v1)
		return trg->v2;
	return NULL;
}

inline msh_edge_t* mtrg_get_opposite_edge(const msh_trg_t *const trg,
					  const msh_vtx_t *const vtx)
{
	if (trg->v1 == vtx)
		return trg->s2;
	if (trg->v2 == vtx)
		return trg->s3;
	if (trg->v3 == vtx)
		return trg->s1;
	return NULL;
}

inline msh_edge_t* mtrg_get_right_edge(const msh_trg_t *const trg, 
				       const msh_vtx_t *const vtx)
{
	if (trg->v1 == vtx)
		return trg->s1;
	if (trg->v2 == vtx)
		return trg->s2;
	if (trg->v3 == vtx)
		return trg->s3;
	return NULL;
}

inline msh_edge_t* mtrg_get_left_edge(const msh_trg_t *const trg, 
				      const msh_vtx_t *const vtx)
{
	if (trg->v1 == vtx)
		return trg->s3;
	if (trg->v2 == vtx)
		return trg->s1;
	if (trg->v3 == vtx)
		return trg->s2;
	return NULL;
}

void mtrg_get_complement_edges(const msh_trg_t *const trg,
			       const msh_edge_t *const edge,
			       msh_edge_t* edge_complement[2])
{
	if (trg->s1 == edge) {
		edge_complement[0] = trg->s2;
		edge_complement[1] = trg->s3;
		return;
	}
	if (trg->s2 == edge) {
		edge_complement[0] = trg->s1;
		edge_complement[1] = trg->s3;
		return;
	}
	if (trg->s3 == edge) {
		edge_complement[0] = trg->s1;
		edge_complement[1] = trg->s2;
		return;
	}
	edge_complement[0] = NULL;
	edge_complement[1] = NULL;
}

inline msh_trg_t* mtrg_get_right_triangle(const msh_trg_t *const trg, 
					  const msh_vtx_t *const vtx)
{
	if (trg->v1 == vtx)
		return trg->t1;
	if (trg->v2 == vtx)
		return trg->t2;
	if (trg->v3 == vtx)
		return trg->t3;
	return NULL;
}

inline msh_trg_t* mtrg_get_left_triangle(const msh_trg_t *const trg, 
					 const msh_vtx_t *const vtx)
{
	if (trg->v1 == vtx)
		return trg->t3;
	if (trg->v2 == vtx)
		return trg->t1;
	if (trg->v3 == vtx)
		return trg->t2;
	return NULL;
}

inline void mtrg_vanish_from_neighbour
            (const msh_trg_t *const restrict trg, 
	     msh_trg_t *const restrict nb_trg)
{
	if (nb_trg->t1 == trg)
		nb_trg->t1 = NULL;
	else if (nb_trg->t2 == trg)
		nb_trg->t2 = NULL;
	else if (nb_trg->t3 == trg)
		nb_trg->t3 = NULL;
}

inline void mtrg_disconnect(const msh_trg_t *const restrict trg)
{
	if (NULL != trg->t1)
		mtrg_vanish_from_neighbour(trg, trg->t1);
	  
	if (NULL != trg->t2)
		mtrg_vanish_from_neighbour(trg, trg->t2);

	if (NULL != trg->t3)
		mtrg_vanish_from_neighbour(trg, trg->t3);
}

inline msh_edge_t* mtrg_get_CCW_edge(const msh_trg_t *const trg,
				     const msh_edge_t *const edge)
{
	if (trg->s1 == edge)
		return trg->s2;
	if (trg->s2 == edge)
		return trg->s3;
	if (trg->s3 == edge)
		return trg->s1;
	return NULL;  
}

inline msh_edge_t* mtrg_get_CW_edge(const msh_trg_t *const trg,
				    const msh_edge_t *const edge)
{
	if (trg->s1 == edge)
		return trg->s3;
	if (trg->s2 == edge)
		return trg->s1;
	if (trg->s3 == edge)
		return trg->s2;
	return NULL;
}

msh_edge_t* medge_get_CW_subsgm(const msh_edge_t *const restrict sgm,
				const msh_vtx_t *const restrict vtx)
{
	msh_edge_t* edge = (msh_edge_t*)sgm;
	do {
		msh_trg_t* restrict trg = (vtx == edge->v1)? edge->t2: edge->t1;
		if (NULL == trg)
			return NULL;
    
		edge = mtrg_get_CCW_edge(trg, edge); /* It is not a mistake */
		if (medge_is_subsgm(edge))
			return edge;
	} while (edge != sgm);
	return (msh_edge_t*)sgm;
}

msh_edge_t* medge_get_CCW_subsgm(const msh_edge_t *const restrict sgm,
				 const msh_vtx_t *const restrict vtx)
{
	msh_edge_t* edge = (msh_edge_t*) sgm;
	do {
		msh_trg_t* restrict trg = (vtx == edge->v1)? edge->t1: edge->t2;
		if (NULL == trg)
			return NULL;
    
		edge = mtrg_get_CW_edge(trg, edge); /* It is not a mistake */
		if (medge_is_subsgm(edge))
			return edge;
	} while (edge != sgm);
	return (msh_edge_t*) sgm;
}

msh_vtx_t *medge_get_partner_vtx(const msh_edge_t *const edge,
				 const msh_vtx_t *const vtx)
{
	msh_vtx_t *out;
	if (edge->v1 == vtx)
		out = edge->v2;
	else if (edge->v2 == vtx)
		out = edge->v1;
	else
		out = NULL;
	return out;
}

void medge_flip_without_dealloc(msh_edge_t* shared_sgm)
/* WARNING: The hash tables of the mesh are not updated. */
{
	/*        v2                   v2
	 *    sA /|\ sD           sA  / \ sD
	 *      / | \                /T1	\
	 *  v3 /T1|  \v4    ===> v3 /_____\ v4
	 *     \  |T2/              \     /
	 *      \ | /                \ T2/
	 *    sB \|/ sC            sB \ / sC
	 *        v1                   v1
	 */
	msh_trg_t* t1 = shared_sgm->t1;
	if (NULL == t1)
		return;

	msh_trg_t* t2 = shared_sgm->t2;
	if (NULL == t2)
		return;

	msh_vtx_t* v1 = shared_sgm->v1;
	msh_vtx_t* v2 = shared_sgm->v2;
	msh_vtx_t* v3 = mtrg_get_opposite_vertex(t1, shared_sgm);
	msh_vtx_t* v4 = mtrg_get_opposite_vertex(t2, shared_sgm);

	msh_edge_t* sA = mtrg_get_left_edge(t1, v3);
	msh_edge_t* sB = mtrg_get_right_edge(t1, v3);
	msh_edge_t* sC = mtrg_get_left_edge(t2, v4);
	msh_edge_t* sD = mtrg_get_right_edge(t2, v4);

	/* Remove and reeinsert shared segment from hash table */
	shared_sgm->v1 = v3;
	shared_sgm->v2 = v4;

	/* Reconnect vertices */
	t1->v1 = v2;
	t1->v2 = v3;
	t1->v3 = v4;

	t2->v1 = v1;
	t2->v2 = v4;
	t2->v3 = v3;

	/* Reconnect segments */
	t1->s1 = sA;
	t1->s2 = shared_sgm;
	t1->s3 = sD;

	t2->s1 = sC;
	t2->s2 = shared_sgm;
	t2->s3 = sB;

	if(sB->t1 == t1)
		sB->t1 = t2;
	else
		sB->t2 = t2;
	if(sD->t1 == t2)
		sD->t1 = t1;
	else
		sD->t2 = t1;

	/* Reconnect triangles */
	t1->t1 = (sA->t1 == t1)?sA->t2:sA->t1;
	t1->t2 = t2;
	t1->t3 = (sD->t1 == t1)?sD->t2:sD->t1;

	t2->t1 = (sC->t1 == t2)?sC->t2:sC->t1;
	t2->t2 = t1;
	t2->t3 = (sB->t1 == t2)?sB->t2:sB->t1;

	/* Update neighbouring triangles */
	if (t2->t3 != NULL) {
		/* Neighbour corresponding to sB */
		if (t2->t3->s1 == sB)
			t2->t3->t1 = t2;
		else if (t2->t3->s2 == sB)
			t2->t3->t2 = t2;
		else if (t2->t3->s3 == sB)
			t2->t3->t3 = t2;
	}
	if (t1->t3 != NULL) {
		/* Neighbour corresponding to sD */
		if (t1->t3->s1 == sD)
			t1->t3->t1 = t1;
		else if (t1->t3->s2 == sD)
			t1->t3->t2 = t1;
		else if (t1->t3->s3 == sD)
			t1->t3->t3 = t1;
	}
}

inline msh_edge_t* mesh_insert_edge(nb_container_t *const ht_edge,
				    const msh_vtx_t *const v1, 
				    const msh_vtx_t *const v2)
{
	msh_edge_t *const restrict new_sgm = calloc(1, sizeof(*new_sgm));
  
	new_sgm->v1 = (msh_vtx_t*)v1;
	new_sgm->v2 = (msh_vtx_t*)v2;
	nb_container_insert(ht_edge, new_sgm);
	return new_sgm;
}
inline msh_edge_t* mesh_exist_edge_guided
                         (nb_container_t *const restrict ht_edge,
			  const msh_vtx_t *const restrict v1,
			  const msh_vtx_t *const restrict v2)
{
	msh_edge_t key_edge;
	key_edge.v1 = (msh_vtx_t*)v1;
	key_edge.v2 = (msh_vtx_t*)v2;
	return nb_container_exist(ht_edge, &key_edge);
}

inline msh_edge_t* mesh_exist_edge(nb_container_t *const restrict ht_edge,
				   const msh_vtx_t *const restrict v1,
				   const msh_vtx_t *const restrict v2)
{
	msh_edge_t *sgm = mesh_exist_edge_guided(ht_edge, v1, v2);
	if (NULL == sgm)
		sgm = mesh_exist_edge_guided(ht_edge, v2, v1);
	return sgm;
}

void mesh_add_triangle(vcn_mesh_t *const mesh, msh_trg_t *const trg)
{
	/* Insert new triangle into the mesh */
	nb_container_insert(mesh->ht_trg, trg);

	/* Connect triangle with the first segment */
	msh_edge_t* sgm = mesh_exist_edge(mesh->ht_edge, trg->v1, trg->v2);
	if (NULL == sgm)
		sgm = mesh_insert_edge(mesh->ht_edge, trg->v1, trg->v2);
	trg->s1 = sgm;
	if (sgm->v1 == trg->v1)
		sgm->t1 = trg;
	else
		sgm->t2 = trg;
	medge_connect_triangles(sgm);

	/* Connect triangle with the second segment */
	sgm = mesh_exist_edge(mesh->ht_edge, trg->v2, trg->v3);
	if (NULL == sgm)
		sgm = mesh_insert_edge(mesh->ht_edge, trg->v2, trg->v3);
	trg->s2 = sgm;
	if (sgm->v1 == trg->v2)
		sgm->t1 = trg;
	else
		sgm->t2 = trg;
	medge_connect_triangles(sgm);

	/* Connect triangle with the third segment */
	sgm = mesh_exist_edge(mesh->ht_edge, trg->v3, trg->v1);
	if (NULL == sgm)
		sgm = mesh_insert_edge(mesh->ht_edge, trg->v3, trg->v1);
	trg->s3 = sgm;
	if (sgm->v1 == trg->v3)
		sgm->t1 = trg;
	else
		sgm->t2 = trg;
	medge_connect_triangles(sgm);
}

void mesh_substract_triangle(vcn_mesh_t *restrict mesh, 
			     msh_trg_t *restrict trg)
{
	/* Remove from hash table */
	nb_container_delete(mesh->ht_trg, trg);

	/* Disconnect from segments */
	msh_edge_t* sgm = trg->s1;
	if (sgm->v1 == trg->v1) 
		sgm->t1 = NULL;
	else
		sgm->t2 = NULL;
	
	if (!medge_is_subsgm(sgm)) {
		if (sgm->t1 == NULL && sgm->t2 == NULL)
			mesh_remove_edge(mesh->ht_edge, sgm->v1, sgm->v2);
	}

	sgm = trg->s2;
	if (sgm->v1 == trg->v2)
		sgm->t1 = NULL;
	else
		sgm->t2 = NULL;

	if (!medge_is_subsgm(sgm)) {
		if (sgm->t1 == NULL && sgm->t2 == NULL)
			mesh_remove_edge(mesh->ht_edge, sgm->v1, sgm->v2);
	}

	sgm = trg->s3;
	if (sgm->v1 == trg->v3)
		sgm->t1 = NULL;
	else
		sgm->t2 = NULL;

	if (!medge_is_subsgm(sgm)) {
		if (sgm->t1 == NULL && sgm->t2 == NULL)
			mesh_remove_edge(mesh->ht_edge, sgm->v1, sgm->v2);
	}

	/* Disconnect from neigbouring triangles */
	mtrg_disconnect(trg);
}

static bool mesh_remove_edge(nb_container_t *const restrict ht_edge,
			     const msh_vtx_t *const restrict v1, 
			     const msh_vtx_t *const restrict v2)
{
	/* Generate Index and Hash Key */
	msh_edge_t key_sgm;
	key_sgm.v1 = (msh_vtx_t*)v1;
	key_sgm.v2 = (msh_vtx_t*)v2;

	key_sgm.attr = NULL;
	msh_edge_t* sgm = nb_container_exist(ht_edge, &key_sgm);
	if (NULL != sgm) {
		nb_container_delete(ht_edge, sgm);
		free(sgm);
		return true;
	}
	return false;
}

inline uint32_t hash_key_edge(const void *const edge_ptr)
{
	const msh_edge_t *const restrict edge = edge_ptr;
	return (uint32_t)
		((int)(edge->v1->x[0] * 73856093) ^ 
		 (int)(edge->v1->x[1] * 19349663) ^
		 (int)(edge->v2->x[0] * 83492791) ^
		 (int)(edge->v2->x[1] * 83492791));
}

inline int8_t compare_edge(const void *const edge1_ptr,
			   const void *const edge2_ptr)
{
	const msh_edge_t *const edge1 = edge1_ptr;
	const msh_edge_t *const edge2 = edge2_ptr;
	return ((edge1->v1 == edge2->v1) && (edge1->v2 == edge2->v2)) ?
	  0:1;
}

msh_trg_t* mesh_locate_vtx(const vcn_mesh_t *const restrict mesh,
			   const msh_vtx_t *const restrict v)
{
	nb_iterator_t *iter = alloca(nb_iterator_get_memsize());
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, mesh->ht_trg);
	msh_trg_t *enveloping_trg = NULL;
	while (nb_iterator_has_more(iter)) {
		const msh_trg_t *trg = nb_iterator_get_next(iter);
		if (vcn_utils2D_pnt_lies_in_trg(trg->v1->x,
						trg->v2->x,
						trg->v3->x,
						v->x)) {
			enveloping_trg = (msh_trg_t*) trg;
			break;
		}
	}
	nb_iterator_finish(iter);
	return enveloping_trg;
}

inline void mesh_get_extern_scale_and_disp(const vcn_mesh_t *const mesh,
					   const double internal[2],
					   double external[2])
{
	external[0] = internal[0] / mesh->scale + mesh->xdisp;
	external[1] = internal[1] / mesh->scale + mesh->ydisp;
}

void mesh_alloc_vtx_ids(vcn_mesh_t * restrict mesh)
{
	vcn_bins2D_iter_t* iter = vcn_bins2D_iter_create();
	vcn_bins2D_iter_set_bins(iter, mesh->ug_vtx);
	int i = 0;
	while (vcn_bins2D_iter_has_more(iter)) {
	  msh_vtx_t* vtx = (msh_vtx_t*) vcn_bins2D_iter_get_next(iter);
		void** attr = malloc(2 * sizeof(*attr));
		uint32_t* id = malloc(sizeof(*id));
		id[0] = i++;
		attr[0] = id;
		attr[1] = vtx->attr;
		vtx->attr = attr;
	}
	vcn_bins2D_iter_destroy(iter);
}

void mesh_free_vtx_ids(vcn_mesh_t *mesh)
{
	vcn_bins2D_iter_t* iter = vcn_bins2D_iter_create();
	vcn_bins2D_iter_set_bins(iter, mesh->ug_vtx);
	while (vcn_bins2D_iter_has_more(iter)) {
		msh_vtx_t* vtx = (msh_vtx_t*) vcn_bins2D_iter_get_next(iter);
		void** attr = vtx->attr;
		vtx->attr = attr[1];
		free(attr[0]);
		free(attr);
	}
	vcn_bins2D_iter_destroy(iter);
}

void mesh_alloc_trg_ids(vcn_mesh_t *mesh)
{
	uint16_t iter_size = nb_iterator_get_memsize();
	nb_iterator_t* trg_iter = alloca(iter_size);
	nb_iterator_init(trg_iter);
	nb_iterator_set_container(trg_iter, mesh->ht_trg);
	uint32_t i = 0;
	while (nb_iterator_has_more(trg_iter)) {
		msh_trg_t* trg = (msh_trg_t*) nb_iterator_get_next(trg_iter);
		char *memblock = malloc(2 * sizeof(void*) + sizeof(uint32_t*));
		void** attr = memblock;
		uint32_t* id = memblock + 2 * sizeof(void*);
		id[0] = i;
		attr[0] = id;
		attr[1] = trg->attr;
		trg->attr = attr;
		i++;
	}
	nb_iterator_finish(trg_iter);
}

void mesh_free_trg_ids(vcn_mesh_t *mesh)
{
	nb_iterator_t* iter = nb_iterator_create();
	nb_iterator_set_container(iter, mesh->ht_trg);
	while (nb_iterator_has_more(iter)) {
		msh_trg_t* trg = (msh_trg_t*) nb_iterator_get_next(iter);
		void** attr = trg->attr;
		trg->attr = attr[1];
		free(attr);
	}
	nb_iterator_destroy(iter);
}


void mesh_split_trg_formed_by_input_sgm(nb_mesh_t *mesh)
{
	/* PENDING */
}
