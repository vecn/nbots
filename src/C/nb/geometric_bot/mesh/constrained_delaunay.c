#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>

#include "nb/container_bot/container.h"
#include "nb/container_bot/iterator.h"
#include "nb/geometric_bot/utils2D.h"
#include "nb/geometric_bot/knn/bins2D.h"
#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/geometric_bot/mesh/dewall.h"

#include "mesh2D_structs.h"

#define POW2(a) ((a)*(a))

static void set_constraints(vcn_mesh_t *mesh,
			    const uint32_t *const input_sgm,
			    uint32_t N_input_sgm);
static void set_new_constraining_sgm(vcn_mesh_t *mesh,
				     const msh_vtx_t *const v1,
				     const msh_vtx_t *const v2,
				     uint32_t sgm_id);
static vcn_container_t* remove_trg_intersecting_sgm(vcn_mesh_t *mesh,
						    const msh_vtx_t *const v1,
						    const msh_vtx_t *const v2);
static void spread_trg_adjacent_to_constrains
                  (vcn_mesh_t *const mesh,
		   const vcn_container_t *const available_vtx,
		   msh_edge_t *const sgm,
		   bool left_side);
static msh_trg_t* create_trg_constrained
                          (vcn_mesh_t *const mesh,
			   const vcn_container_t *const vertices,
			   msh_edge_t *const base_sgm,
			   const msh_vtx_t *const base_v1);

static bool trg_is_constrained_exahustive_search
                    (msh_edge_t** input_sgm, uint32_t N_input_sgm,
		     const msh_vtx_t *const trg_v1, 
		     const msh_vtx_t *const trg_v2,
		     const msh_vtx_t *const trg_v3,
		     const vcn_container_t *const vertices);
static bool vtx_is_constrained
                    (msh_edge_t** input_sgm, uint32_t N_input_sgm,
		     const msh_vtx_t *const trg_v1,
		     const msh_vtx_t *const trg_v2,
		     const msh_vtx_t *const trg_v3,
		     const msh_vtx_t *const encroaching_vtx);
static bool trg_is_constrained(const vcn_mesh_t *const mesh,
			       const msh_vtx_t *const trg_v1, 
			       const msh_vtx_t *const trg_v2,
			       const msh_vtx_t *const trg_v3,
			       const vcn_container_t *const vertices);
static void link_constraining_sgm(msh_edge_t *sgm);
static bool are_intersecting_edges(const vcn_mesh_t *const mesh,
				   const msh_vtx_t *const v1, 
				   const msh_vtx_t *const v2);

void vcn_mesh_get_constrained_delaunay(vcn_mesh_t *mesh,
				       uint32_t N_vertices,
				       const double *const vertices,
				       uint32_t N_segments,
				       const uint32_t *const segments)
{
	vcn_mesh_get_delaunay(mesh, N_vertices, vertices);
	set_constraints(mesh, segments, N_segments);
}

static void set_constraints(vcn_mesh_t *mesh,
			    const uint32_t *const input_sgm,
			    uint32_t N_input_sgm)
{
	if (0 < N_input_sgm) {
		mesh->N_input_sgm = N_input_sgm;
		mesh->input_sgm = 
			calloc(N_input_sgm, sizeof(*(mesh->input_sgm)));
		for (uint32_t i = 0; i < N_input_sgm; i++) {
			msh_vtx_t* v1 = mesh->input_vtx[input_sgm[i * 2]];
			msh_vtx_t* v2 = mesh->input_vtx[input_sgm[i*2+1]];
			msh_edge_t* sgm = 
				mesh_exist_edge(mesh->ht_edge, v1, v2);
			if (NULL == sgm) {
				set_new_constraining_sgm(mesh, v1, v2, i);
			} else {
				link_constraining_sgm(sgm);
				mesh->input_sgm[i] = sgm;
				medge_set_as_subsgm(sgm, i, NULL, NULL);
			}
		}
	}
}

static void set_new_constraining_sgm(vcn_mesh_t *mesh,
				     const msh_vtx_t *const v1,
				     const msh_vtx_t *const v2,
				     uint32_t sgm_id)
{
	vcn_container_t *vertices = remove_trg_intersecting_sgm(mesh, v1, v2);

	if (vcn_container_is_not_empty(vertices)) {
		msh_edge_t *sgm = mesh_insert_edge(mesh->ht_edge, v1, v2);
		mesh->input_sgm[sgm_id] = sgm;
		medge_set_as_subsgm(sgm, sgm_id, NULL, NULL);
		/* Triangulate left side */
		spread_trg_adjacent_to_constrains(mesh, vertices, sgm, true);
		/* Triangulate right side */
		spread_trg_adjacent_to_constrains(mesh, vertices, sgm, false);
	}
	vcn_container_destroy(vertices);
}

static vcn_container_t* remove_trg_intersecting_sgm(vcn_mesh_t *mesh,
						    const msh_vtx_t *const v1,
						    const msh_vtx_t *const v2)
{
	vcn_container_t *intersected_trg = 
		vcn_container_create(NB_CONTAINER_QUEUE);
	vcn_iterator_t *iter = vcn_iterator_create();
	vcn_iterator_set_container(iter, mesh->ht_trg);
	while (vcn_iterator_has_more(iter)) {
		const msh_trg_t* trg = vcn_iterator_get_next(iter);
		if (vcn_utils2D_sgm_intersects_trg(trg->v1->x,
						   trg->v2->x,
						   trg->v3->x,
						   v1->x, v2->x)) {
			vcn_container_insert(intersected_trg, trg);
		}
	}
	vcn_iterator_destroy(iter);

	vcn_container_t *vertices = 
		vcn_container_create(NB_CONTAINER_SORTED);
	while (vcn_container_is_not_empty(intersected_trg)) {
		msh_trg_t *trg = vcn_container_delete_first(intersected_trg);
		vcn_container_insert(vertices, trg->v1);
		vcn_container_insert(vertices, trg->v2);
		vcn_container_insert(vertices, trg->v3);
		mesh_substract_triangle(mesh, trg);
		free(trg);
	}
	vcn_container_destroy(intersected_trg);
	return vertices;
}

static void spread_trg_adjacent_to_constrains
                  (vcn_mesh_t *const restrict mesh,
		   const vcn_container_t *const restrict available_vtx,
		   msh_edge_t *const restrict sgm,
		   bool left_side)
{
	msh_vtx_t* restrict v1 = (left_side)?sgm->v1:sgm->v2;
	msh_vtx_t* restrict v2 = (left_side)?sgm->v2:sgm->v1;

	msh_trg_t* new_trg =
		create_trg_constrained(mesh, available_vtx, sgm, v1);
	if (NULL != new_trg) {
		mesh->do_after_insert_trg(mesh);

  		/* Process Counter Clock Wise segment */
		msh_edge_t* nb_sgm = mtrg_get_CCW_edge(new_trg, sgm);
	
		bool process_left_side = true;
		if (nb_sgm->v1 == v2)
			process_left_side = false;

		spread_trg_adjacent_to_constrains(mesh, available_vtx, 
						  nb_sgm, process_left_side);

		/* Process Clock Wise segment */
		nb_sgm = mtrg_get_CW_edge(new_trg, sgm);

		process_left_side = true;
		if (nb_sgm->v2 == v1)
			process_left_side = false;

		spread_trg_adjacent_to_constrains(mesh, available_vtx, 
						  nb_sgm, process_left_side);
	}

}


static msh_trg_t* create_trg_constrained
                          (vcn_mesh_t *const restrict mesh,
			   const vcn_container_t *const restrict vertices,
			   msh_edge_t *const restrict base_sgm,
			   const msh_vtx_t *const restrict base_v1)
{
	msh_edge_t* sgm1 = NULL; /* Dummy initialization */
	msh_edge_t* sgm2 = NULL; /* Dummy initialization */
	const msh_vtx_t *const restrict v1 = base_v1;
	const msh_vtx_t *const restrict v2 = 
		(base_sgm->v1 == base_v1)?base_sgm->v2:base_sgm->v1;

	msh_vtx_t* restrict vtx_near = NULL;
	bool vtx_found = false;
	double min_circumradius = 0; /* Must be initialized with zero */
	vcn_iterator_t* iter = vcn_iterator_create();
	vcn_iterator_set_container(iter, vertices);
	while (vcn_iterator_has_more(iter)) {
		const msh_vtx_t *const v3 =
			(msh_vtx_t*)vcn_iterator_get_next(iter);

		/* Check if it is not the same vertex */
		if (v3 == v1 || v3 == v2)
			continue;
    
		/* Check for positive area */
		const double Sk = vcn_utils2D_get_2x_trg_area(v1->x, v2->x, v3->x);
		if (Sk < NB_GEOMETRIC_TOL)
			continue;

		/* Check if the segments already have a triangle */
		msh_edge_t *const sgm1_i = mesh_exist_edge(mesh->ht_edge, v2, v3);
		if (medge_has_a_triangle_to_the_left(sgm1_i, v2))
			continue;
		msh_edge_t *const sgm2_i = mesh_exist_edge(mesh->ht_edge, v3, v1);
		if (medge_has_a_triangle_to_the_left(sgm2_i, v3))
			continue;

		/* Compute distance between vertices */
		const double L1 = vcn_utils2D_get_dist(v1->x, v2->x);
		const double L2 = vcn_utils2D_get_dist(v3->x, v2->x);
		const double L3 = vcn_utils2D_get_dist(v3->x, v1->x);
    
		/* Compute circumradius and circumcenter*/
		const double circumradius = (L1*L2*L3)/(2.0*Sk);
    
		/* Select the vertex with the minimum circumradius */
		if (!vtx_found || min_circumradius > circumradius-NB_GEOMETRIC_TOL) {
			/* Check if the triangle accomplish the Delaunay condition */
			if (!trg_is_constrained(mesh,
									    v1, v2, v3,
									    vertices))
				continue;
			/* Check if there are not intersections with segments */
			if (min_circumradius < circumradius + NB_GEOMETRIC_TOL) {
				if (are_intersecting_edges(mesh, v1, v3))
					continue;
				if (are_intersecting_edges(mesh, v2, v3))
					continue;
			}
			/* Set first vlues */
			min_circumradius = circumradius;
			vtx_near = (msh_vtx_t*)v3; /* Casting from const to non-const pointer
						    * in order to assign it to the new triangle,
						    * but with the compromise to do not modify
						    * its value */
			/* Update segments and its status */
			sgm1 = sgm1_i;
			sgm2 = sgm2_i;
			vtx_found = true;
		}
	}
	vcn_iterator_destroy(iter);

	if (!vtx_found)
		return NULL;
	/* Insert vertex into the mesh */
	const msh_vtx_t *const restrict v3 = vtx_near;

	/* Create and insert triangle */
	msh_trg_t *const restrict new_trg = calloc(1, sizeof(*new_trg));

	new_trg->v1 = (msh_vtx_t*)v1; /* Casting from const to non-const pointer */
	new_trg->v2 = (msh_vtx_t*)v2; /* in order to assign it to the new        */
	new_trg->v3 = (msh_vtx_t*)v3; /* triangle, but with the compromise to do
				       * not modify its value */

	vcn_container_insert(mesh->ht_trg, new_trg);
	/* Connect triangle with the first segment */
	medge_connect_triangle(base_sgm, new_trg, v1, v2);
	/* Connect triangle with the second segment */
	if (NULL == sgm1)
		sgm1 = mesh_insert_edge(mesh->ht_edge, v2, v3);
	medge_connect_triangle(sgm1, new_trg, v2, v3);
	/* Connect triangle with the third segment */
	if (NULL == sgm2)
		sgm2 = mesh_insert_edge(mesh->ht_edge, v3, v1);
	medge_connect_triangle(sgm2, new_trg, v3, v1);
	return new_trg;
}

static bool trg_is_constrained_exahustive_search
                    (msh_edge_t** input_sgm, uint32_t N_input_sgm,
		     const msh_vtx_t *const restrict trg_v1, 
		     const msh_vtx_t *const restrict trg_v2,
		     const msh_vtx_t *const restrict trg_v3,
		     const vcn_container_t *const restrict vertices)
/* List Naive Search */
{
	/* Check the Delaunay condition */
	vcn_iterator_t *const restrict iter = vcn_iterator_create();
	vcn_iterator_set_container(iter, vertices);
	while (vcn_iterator_has_more(iter)) {
		const msh_vtx_t *const restrict vtx = 
			vcn_iterator_get_next(iter);
		if (vtx == trg_v1 || vtx == trg_v2 || vtx == trg_v3)
			continue;
		if (vcn_utils2D_pnt_lies_strictly_in_circumcircle(trg_v1->x,
								  trg_v2->x,
								  trg_v3->x,
								  vtx->x)) {
			/* Check for constraints */
			if (!vtx_is_constrained(input_sgm, N_input_sgm,
						trg_v1, trg_v2, trg_v3, vtx)) {
				vcn_iterator_destroy(iter);
				return false;
			}
		}
	}
	vcn_iterator_destroy(iter);
	return true;
}

static bool vtx_is_constrained
                    (msh_edge_t** input_sgm, uint32_t N_input_sgm,
		     const msh_vtx_t *const restrict trg_v1,
		     const msh_vtx_t *const restrict trg_v2,
		     const msh_vtx_t *const restrict trg_v3,
		     const msh_vtx_t *const restrict encroaching_vtx)
/* Naive search */
{
	/* Check for constraints */
	for (uint32_t i=0; i < N_input_sgm; i++) {
		msh_edge_t* sgm = input_sgm[i];
		while (NULL != sgm) {
			if(vcn_utils2D_are_sgm_intersected(trg_v1->x,
							   encroaching_vtx->x,
							   sgm->v1->x,
							   sgm->v2->x,
							   NULL, NULL))
				return true;
			if(vcn_utils2D_are_sgm_intersected(trg_v2->x,
							   encroaching_vtx->x,
							   sgm->v1->x,
							   sgm->v2->x,
							   NULL, NULL))
				return true;
			if(vcn_utils2D_are_sgm_intersected(trg_v3->x,
							   encroaching_vtx->x,
							   sgm->v1->x,
							   sgm->v2->x,
							   NULL, NULL))
				return true;
			sgm = medge_subsgm_next(sgm);
		}
	}
	return false;
}

static bool trg_is_constrained(const vcn_mesh_t *const restrict mesh,
			       const msh_vtx_t *const restrict trg_v1, 
			       const msh_vtx_t *const restrict trg_v2,
			       const msh_vtx_t *const restrict trg_v3,
			       const vcn_container_t *const restrict vertices)
{
	if (vcn_container_get_length(vertices) < 100)
		return trg_is_constrained_exahustive_search
			(mesh->input_sgm,
			 mesh->N_input_sgm, 
			 trg_v1, trg_v2, trg_v3, vertices);

	const double circumradius =
	  vcn_utils2D_get_circumradius(trg_v1->x, trg_v2->x, trg_v3->x);
	const int layer = 1 +
		(int)(circumradius/vcn_bins2D_get_size_of_bins(mesh->ug_vtx));
	if (POW2(2 * layer + 1) >=  vcn_bins2D_get_N_bins(mesh->ug_vtx))
		return trg_is_constrained_exahustive_search(mesh->input_sgm,
							    mesh->N_input_sgm,
							    trg_v1, trg_v2,
							    trg_v3, vertices);

	/* Get vertices encroaching the circumcircle */
	double circumcenter[2];
	vcn_utils2D_get_circumcenter(trg_v1->x, trg_v2->x, trg_v3->x,
				     circumcenter);
	vcn_container_t *const restrict l_inside_vtx = 
		vcn_bins2D_get_points_inside_circle(mesh->ug_vtx,
						    circumcenter,
						    circumradius);
	bool is_Constrained_Delaunay = true;
	if (vcn_container_is_not_empty(l_inside_vtx)) {
		is_Constrained_Delaunay = 
			trg_is_constrained_exahustive_search(mesh->input_sgm, 
							     mesh->N_input_sgm, 
							     trg_v1, trg_v2,
							     trg_v3,
							     l_inside_vtx);
	}
	vcn_container_destroy(l_inside_vtx);
	return is_Constrained_Delaunay;
}

static void link_constraining_sgm(msh_edge_t *sgm)
{
	msh_trg_t* t1 = sgm->t1;
	msh_trg_t* t2 = sgm->t2;
	msh_vtx_t *v1 = sgm->v1;
	if(t1 != NULL){
		if(v1 == t1->v1)
			t1->s1 = sgm;
		else if(v1 == t1->v2)
			t1->s2 = sgm;
		else if(v1 == t1->v3)
			t1->s3 = sgm;
	}
	if(t2 != NULL){
		if(v1 == t2->v2)
			t2->s1 = sgm;
		else if(v1 == t2->v3)
			t2->s2 = sgm;
		else if(v1 == t2->v1)
			t2->s3 = sgm;
	}
}

static bool are_intersecting_edges(const vcn_mesh_t *const restrict mesh,
				   const msh_vtx_t *const restrict v1, 
				   const msh_vtx_t *const restrict v2)
/* Exahustive search */
{
	vcn_iterator_t *const restrict iter = vcn_iterator_create();
	vcn_iterator_set_container(iter, mesh->ht_edge);

	while (vcn_iterator_has_more(iter)) {
		const msh_edge_t *const restrict edge =
			(msh_edge_t*)vcn_iterator_get_next(iter);

		if (v1 == edge->v1 || v1 == edge->v2 ||
		    v2 == edge->v1 || v2 == edge->v2)
			continue;
		if (vcn_utils2D_are_sgm_intersected(edge->v1->x,
						    edge->v2->x,
						    v1->x,
						    v2->x,
						    NULL, NULL)) {
			vcn_iterator_destroy(iter);
			return true;
		}
	}
	vcn_iterator_destroy(iter);
	return false;
}
