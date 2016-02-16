#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "nb/container_bot.h"
#include "nb/geometric_bot/point2D.h"
#include "nb/geometric_bot/utils2D.h"
#include "nb/geometric_bot/knn/bins2D.h"
#include "nb/geometric_bot/knn/bins2D_iterator.h"
#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/geometric_bot/mesh/dewall.h"
#include "nb/geometric_bot/mesh/modules2D/graph_generator.h"
#include "nb/geometric_bot/mesh/elements2D/polygons.h"

#include "../mesh2D_structs.h"

#define POW2(a) ((a)*(a))

typedef struct {
	msh_edge_t* sgm;
	msh_vtx_t* v1;
	uint32_t N_vertices; /* 1 or 2 vertices */
	vcn_point2D_t* vertices[2];
	bool vtx_inserted[2];
} voronoi_boundary_sgm_t;

static vcn_point2D_t* voronoi_insert_global(vcn_bins2D_t* grid, vcn_point2D_t* point);

static void voronoi_get_patched_circumcenter(double* v1,
					     double* v2,
					     double* v3,
					     /* Output */
					     double* circumcenter);
static nb_container_t* voronoi_get_boundary_segments_with_vertices
                                        (const vcn_mesh_t* const mesh);
		
static int8_t voronoi_vtx_compare(const void* const A,
				  const void *const B);

static inline int compare_voronoi_sgm(const void *const sgmA, 
				      const void *const sgmB);

static void mesh_Lloyd_iteration(msh_vtx_t** vertices, uint32_t N_vertices,
				 msh_edge_t** segments, uint32_t N_segments,
				 msh_trg_t** triangles, uint32_t N_triangles,
				 /* NULL for a constant density */
				 double (*density)(double*),
				 const nb_container_t *const avl_boundary_sgm,
				 uint32_t max_iter);

vcn_mshpoly_t* vcn_mesh_get_mshpoly(const vcn_mesh_t *const mesh,
				    bool include_adjacencies,
				    bool central_voronoi,
				    uint32_t central_voronoi_max_iter,
				    /* NULL for a constant density */
				    double (*central_voronoi_density)(double*),
				    uint32_t* (*labeling)
				    		(const vcn_graph_t *const))
{
	/* Allocate Voronoi mesh */
	vcn_mshpoly_t* voronoi = calloc(1, sizeof(*voronoi));

	/* Clone the mesh if the central Voronoi is selected
	 * to allow the vertices modification */
	vcn_mesh_t* aux_mesh = (vcn_mesh_t*) mesh;
	if(central_voronoi)
		aux_mesh = vcn_mesh_clone(mesh);

	/* Put vertices, segments and triangles into
	 * arrays, because the hash key generated is
	 * based on the position of the vertices, which
	 * is changing,  and they can not be accessed 
	 * via hash table. 
	 */
	uint32_t N_vertices = vcn_bins2D_get_length(aux_mesh->ug_vtx);
	msh_vtx_t** vertices =  malloc(2 * voronoi->N_polygons *
				       sizeof(*vertices));
	uint32_t N_segments = nb_container_get_length(aux_mesh->ht_edge);
	msh_edge_t** segments = malloc(N_segments * sizeof(*segments));
	uint32_t N_triangles = nb_container_get_length(aux_mesh->ht_trg);
	msh_trg_t** triangles = malloc(N_triangles * sizeof(*triangles));

	vcn_bins2D_iter_t* iter = vcn_bins2D_iter_create();
	vcn_bins2D_iter_set_bins(iter, aux_mesh->ug_vtx);
	int i = 0;
	while (vcn_bins2D_iter_has_more(iter)) {
		msh_vtx_t* vtx = (msh_vtx_t*) vcn_bins2D_iter_get_next(iter);
		/* Allocate the ID */
		void **attr = malloc(2 * sizeof(*attr));
		uint32_t *id = malloc(sizeof(*id));
		id[0] = i++;
		attr[0] = id;
		attr[1] = vtx->attr;
		vtx->attr = attr;
	}
	nb_iterator_t *sgm_iter = nb_iterator_create();
	nb_iterator_set_container(sgm_iter, aux_mesh->ht_edge);
	i = 0;
	while (nb_iterator_has_more(sgm_iter)) {
		msh_edge_t* sgm = (msh_edge_t*)nb_iterator_get_next(sgm_iter);
		segments[i++] = sgm;
	}
	nb_iterator_destroy(sgm_iter);

	nb_iterator_t *trg_iter = nb_iterator_create();
	nb_iterator_set_container(trg_iter, aux_mesh->ht_trg);
	i = 0;
	while (nb_iterator_has_more(trg_iter)) {
		msh_trg_t* trg = (msh_trg_t*)nb_iterator_get_next(trg_iter);
		/* Allocate the ID */
		void **attr = malloc(2 * sizeof(*attr));
		uint32_t *id = malloc(sizeof(*id));
		id[0] = i;
		attr[0] = id;
		attr[1] = trg->attr;
		trg->attr = attr;
		/* Put the triangle into the array */
		triangles[i++] = trg;
	}
	nb_iterator_destroy(trg_iter);

	uint32_t *perm = NULL; /* Dummy initialization */
	if (NULL != labeling) {
		vcn_graph_t *graph = vcn_mesh_create_vtx_graph(mesh);
		perm = labeling(graph);
		vcn_graph_destroy(graph);
	}

	/* Allocate Voronoi data */
	voronoi->N_polygons = N_vertices;
	voronoi->centroids = 
		malloc(2 * voronoi->N_polygons * sizeof(*(voronoi->centroids)));
	voronoi->N_vertices_forming_polygons =
		calloc(voronoi->N_polygons,
		       sizeof(*(voronoi->N_vertices_forming_polygons)));
	voronoi->vertices_forming_polygons =
		malloc(voronoi->N_polygons * 
		       sizeof(*(voronoi->vertices_forming_polygons)));

	/* Set the vertices into the array */
	vcn_bins2D_iter_restart(iter);
	while (vcn_bins2D_iter_has_more(iter)) {
		msh_vtx_t* vtx = (msh_vtx_t*)vcn_bins2D_iter_get_next(iter);
		uint32_t* id = (uint32_t*)((void**)vtx->attr)[0];
		if (NULL != labeling)
			id[0] = perm[id[0]];
		/* Put the vertex in the array */
		vertices[id[0]] = vtx;
	}
	vcn_bins2D_iter_destroy(iter);

	/* Free permutation vector used for labeling */
	if (NULL != labeling)
		free(perm);

	/* Compute vertices of polygons which lies on the boundary */
	nb_container_t* avl_boundary_sgm = 
		voronoi_get_boundary_segments_with_vertices(aux_mesh);

	/* Lloyd Iteration */
	if(central_voronoi)
		mesh_Lloyd_iteration(vertices, N_vertices,
				     segments, N_segments,
				     triangles, N_triangles,
				     central_voronoi_density,
				     avl_boundary_sgm, 
				     central_voronoi_max_iter);

	/* Array to access the triangles connected to each vertex */
	msh_trg_t** trg_x_vtx = calloc(N_vertices, sizeof(*trg_x_vtx));

	/* Array with the nodes conforming each Voronoi polygon */
	nb_container_t** node_x_vtx = malloc(N_vertices * sizeof(*node_x_vtx));
	for (uint32_t i = 0; i < N_vertices; i++)
		node_x_vtx[i] = nb_container_create(NB_STACK);

	/* Array to store the circumcenter of the triangles */
	vcn_point2D_t** trg_center = calloc(N_triangles, sizeof(*trg_center));

	/* Relate the triangle with its vertices */  
	for (uint32_t i = 0; i < N_triangles; i++) {
		msh_trg_t* trg = triangles[i];
		trg_x_vtx[((uint32_t*)((void**)trg->v1->attr)[0])[0]] = trg;
		trg_x_vtx[((uint32_t*)((void**)trg->v2->attr)[0])[0]] = trg;
		trg_x_vtx[((uint32_t*)((void**)trg->v3->attr)[0])[0]] = trg;
	}


	/* Compute vertices of polygons */
	vcn_bins2D_t* voronoi_vtx = vcn_bins2D_create(1.0);
	for (uint32_t i = 0; i < voronoi->N_polygons; i++) {
		msh_trg_t* first_trg = trg_x_vtx[i];
		/* Check if the flower is complete */
		msh_trg_t* trg = first_trg;
		while (NULL != trg) {
			/* Next triangle */
			trg = mtrg_get_left_triangle(trg, vertices[i]);
			if(trg == first_trg)
				break;      
		}
		if (trg != first_trg) {
			trg = first_trg;
			while (NULL != trg) {
				first_trg = trg;
				/* Prev triangle */
				trg = mtrg_get_right_triangle(trg, vertices[i]);
			}
			msh_edge_t* ext_sgm = first_trg->s1;
			if (first_trg->v2 == vertices[i])
				ext_sgm = first_trg->s2;
			else if (first_trg->v3 == vertices[i])
				ext_sgm = first_trg->s3;
			/* Fix boundary */
			if (medge_is_subsgm(ext_sgm)) {
				/* Insert boundary vertex */
				voronoi_boundary_sgm_t aux_bsgm;
				aux_bsgm.sgm = ext_sgm;
				voronoi_boundary_sgm_t* bsgm =
					nb_container_exist(avl_boundary_sgm, &aux_bsgm);
      
				/* Insert centroid */
				vcn_point2D_t *centroid = bsgm->v1;
				vcn_point2D_t* centroid_cpy = centroid;
				centroid = voronoi_insert_global(voronoi_vtx, centroid_cpy);
				if (centroid != centroid_cpy)
					vcn_point2D_destroy(centroid_cpy);
				nb_container_insert(node_x_vtx[i], centroid);
				if (bsgm->N_vertices > 0) {
					/* Case 1 & 2 */
					if (!bsgm->vtx_inserted[0]) {
						vcn_bins2D_insert(voronoi_vtx, bsgm->vertices[0]);
						bsgm->vtx_inserted[0] = true;
					}
					nb_container_insert(node_x_vtx[i], bsgm->vertices[0]);
					if (bsgm->N_vertices == 2) {
						/* Case 2 */
						/* Next triangle */
						msh_trg_t* next_trg = first_trg->t3;
						if (first_trg->v2 == vertices[i])
							next_trg = first_trg->t1;
						else if (first_trg->v3 == vertices[i])
							next_trg = first_trg->t2;
						if (next_trg == NULL) {
							/* Corner triangle, Case 1, next boundary */
							msh_edge_t* corner_ext_sgm = first_trg->s3;
							if (first_trg->v2 == vertices[i])
								corner_ext_sgm = first_trg->s1;
							else if (first_trg->v3 == vertices[i])
								corner_ext_sgm = first_trg->s2;
							aux_bsgm.sgm = corner_ext_sgm;
							bsgm = nb_container_exist(avl_boundary_sgm, &aux_bsgm);
							if (!bsgm->vtx_inserted[0]) {
								vcn_bins2D_insert(voronoi_vtx, bsgm->vertices[0]);
								bsgm->vtx_inserted[0] = true;
							}
							nb_container_insert(node_x_vtx[i], bsgm->vertices[0]);
							continue; /* Iteration over Voronoi polygons */
						} else {
							first_trg = next_trg;
						}
					}
				}
			}
		}
		/* Start composition of polygon */
		trg = first_trg;
		msh_trg_t* last_trg = NULL;
		
		while (NULL != trg) {
			uint32_t* trg_id = (uint32_t*)((void**)trg->attr)[0];
			msh_edge_t* op_sgm = mtrg_get_opposite_edge(trg, vertices[i]);
			voronoi_boundary_sgm_t* bsgm = NULL; /* Dummy initialization */
			bool insert_circumcenter_flag = false;
			if (medge_is_subsgm(op_sgm)) {
				if (op_sgm->t1 == NULL || op_sgm->t2 == NULL) {
					/* Insert boundary vertex */
					voronoi_boundary_sgm_t aux_bsgm;
					aux_bsgm.sgm = op_sgm;
					bsgm =
						nb_container_exist(avl_boundary_sgm, &aux_bsgm);
					if (bsgm->N_vertices == 0 || bsgm->N_vertices == 1)
						insert_circumcenter_flag = true;
				} else {
					insert_circumcenter_flag = true;
				}
			} else {
				insert_circumcenter_flag = true;
			}
      
			if (insert_circumcenter_flag) {
				if (NULL == trg_center[trg_id[0]]) {
					/* Compute circumcenter */
					vcn_point2D_t* circumcenter = vcn_point2D_create();
					voronoi_get_patched_circumcenter(trg->v1->x,
									 trg->v2->x,
									 trg->v3->x,
									 circumcenter->x);
					vcn_point2D_t* circumcenter_cpy = circumcenter;
					circumcenter = voronoi_insert_global(voronoi_vtx, circumcenter_cpy);
					if (circumcenter != circumcenter_cpy)
						vcn_point2D_destroy(circumcenter_cpy);
					trg_center[trg_id[0]] = circumcenter;
				}
				if(!nb_container_exist(node_x_vtx[i], trg_center[trg_id[0]]))
					nb_container_insert(node_x_vtx[i], trg_center[trg_id[0]]);
			} else {
				if (!bsgm->vtx_inserted[0]) {
					vcn_bins2D_insert(voronoi_vtx, bsgm->vertices[0]);
					bsgm->vtx_inserted[0] = true;
				}
				if (!bsgm->vtx_inserted[1]) {
					vcn_bins2D_insert(voronoi_vtx, bsgm->vertices[1]);
					bsgm->vtx_inserted[1] = true;
				}
				nb_container_insert(node_x_vtx[i], bsgm->vertices[0]);
				nb_container_insert(node_x_vtx[i], bsgm->vertices[1]);
			}

			/* Next triangle */
			last_trg = trg;
			trg = mtrg_get_left_triangle(trg, vertices[i]);
			if (trg == first_trg)
				break;      
		}

		if (NULL == trg) {
			/* Fix boundary */
			msh_edge_t* ext_sgm = last_trg->s3;
			if (last_trg->v2 == vertices[i])
				ext_sgm = last_trg->s1;
			else if (last_trg->v3 == vertices[i])
				ext_sgm = last_trg->s2;
			if (medge_is_subsgm(ext_sgm)) {
				/* Insert boundary vertex */
				voronoi_boundary_sgm_t aux_bsgm;
				aux_bsgm.sgm = ext_sgm;
				voronoi_boundary_sgm_t* bsgm =
					nb_container_exist(avl_boundary_sgm, &aux_bsgm);

				if (bsgm->N_vertices == 1) {
					/* Case 1 */
					if (!bsgm->vtx_inserted[0]) {
						vcn_bins2D_insert(voronoi_vtx, bsgm->vertices[0]);
						bsgm->vtx_inserted[0] = true;
					}
					nb_container_insert(node_x_vtx[i], bsgm->vertices[0]);
				} else if (bsgm->N_vertices == 2) {
					/* Case 2 */
					uint32_t* last_trg_id = (uint32_t*)((void**)last_trg->attr)[0];
					vcn_point2D_t* last_vtx = nb_container_delete_first(node_x_vtx[i]);
					vcn_bins2D_delete(voronoi_vtx, last_vtx);
					free(trg_center[last_trg_id[0]]);
					trg_center[last_trg_id[0]] = NULL;

					if (!bsgm->vtx_inserted[1]) {
						vcn_bins2D_insert(voronoi_vtx, bsgm->vertices[1]);
						bsgm->vtx_inserted[1] = true;
					}
					nb_container_insert(node_x_vtx[i], bsgm->vertices[1]);
				}
			}
		}
	}
	free(trg_center);

	/* Store local and global info of vertices and centroids*/
	voronoi->N_vertices = vcn_bins2D_get_length(voronoi_vtx);
	voronoi->vertices = malloc(2 * voronoi->N_vertices * sizeof(*(voronoi->vertices)));
	uint32_t vtx_counter = 0;
	for (uint32_t i=0; i < voronoi->N_polygons; i++) {
		/* Store centroids */
		voronoi->centroids[i * 2] = 
			vertices[i]->x[0]/aux_mesh->scale + aux_mesh->xdisp;
		voronoi->centroids[i*2+1] = 
			vertices[i]->x[1]/aux_mesh->scale + aux_mesh->ydisp;
		/* Store vertices (boundaries of polygons) */
		voronoi->N_vertices_forming_polygons[i] = nb_container_get_length(node_x_vtx[i]);
		voronoi->vertices_forming_polygons[i] =
			malloc(voronoi->N_vertices_forming_polygons[i] *
			       sizeof(*(voronoi->vertices_forming_polygons[i])));
		for (uint32_t j = 0; j < voronoi->N_vertices_forming_polygons[i]; j++) {
			vcn_point2D_t* vtx = nb_container_delete_first(node_x_vtx[i]);
			uint32_t* id;
			if (NULL == vtx->attr) {
				id = malloc(sizeof(*id));
				id[0] = vtx_counter++;
				vtx->attr = id;

				voronoi->vertices[id[0] * 2] = 
					vtx->x[0]/aux_mesh->scale + aux_mesh->xdisp;
				voronoi->vertices[id[0]*2+1] = 
					vtx->x[1]/aux_mesh->scale + aux_mesh->ydisp;
			} else {
				id = (uint32_t*) vtx->attr;
			}
			voronoi->vertices_forming_polygons[i][j] = id[0];
		}
		/* Free memory */
		nb_container_destroy(node_x_vtx[i]);
	}

	/* Get adjacencies */
	if (include_adjacencies) {
		/* Create adjacencies matrix */
		voronoi->N_adjacencies = 
			calloc(N_vertices, sizeof(*(voronoi->N_adjacencies)));
		voronoi->adjacencies = 
			malloc(N_vertices * sizeof(*(voronoi->adjacencies)));
  
		for (uint32_t i = 0; i < N_segments; i++) {
			msh_edge_t* sgm = segments[i];
			uint32_t idx1 = ((uint32_t*)((void**)sgm->v1->attr)[0])[0];
			uint32_t idx2 = ((uint32_t*)((void**)sgm->v2->attr)[0])[0];
			voronoi->N_adjacencies[idx1] += 1;
			voronoi->N_adjacencies[idx2] += 1;
		}

		for (uint32_t i = 0; i < N_vertices; i++)
			voronoi->adjacencies[i] =
				malloc(voronoi->N_adjacencies[i] * 
				       sizeof(*(voronoi->adjacencies[i])));

		uint32_t* adj_matrix_next_idx = calloc(N_vertices,
						       sizeof(*adj_matrix_next_idx));

		for (uint32_t i = 0; i < N_segments; i++) {
			msh_edge_t* sgm = segments[i];
			uint32_t idx1 = ((uint32_t*)((void**)sgm->v1->attr)[0])[0];
			uint32_t idx2 = ((uint32_t*)((void**)sgm->v2->attr)[0])[0];
			voronoi->adjacencies[idx1][adj_matrix_next_idx[idx1]] = idx2;
			voronoi->adjacencies[idx2][adj_matrix_next_idx[idx2]] = idx1;
			adj_matrix_next_idx[idx1] += 1;
			adj_matrix_next_idx[idx2] += 1;
		}
		free(adj_matrix_next_idx);
	}

	/* Free memory */
	free(node_x_vtx);
	free(trg_x_vtx);
	vcn_bins2D_set_attribute_destroyer(voronoi_vtx, free);
	vcn_bins2D_destroy(voronoi_vtx);
	nb_container_set_destroyer(avl_boundary_sgm, free);
	nb_container_destroy(avl_boundary_sgm);

	for (uint32_t i = 0; i < N_vertices; i++) {
		msh_vtx_t* vtx = vertices[i];
		void** attr = (void**)vtx->attr;
		vtx->attr = attr[1];
		free(attr[0]);
		free(attr);
	}

	for (uint32_t i = 0; i < N_triangles; i++) {
		msh_trg_t* trg = triangles[i];
		void** attr = (void**)trg->attr;
		trg->attr = attr[1];
		free(attr[0]);
		free(attr);
	}

	free(vertices);
	free(segments);
	free(triangles);

	if (central_voronoi)
		vcn_mesh_destroy(aux_mesh);

	/* Return Voronoi mesh */
	return voronoi;
}

void vcn_mshpoly_destroy(vcn_mshpoly_t* voronoi)
{
	free(voronoi->centroids);
	free(voronoi->vertices);
	free(voronoi->N_vertices_forming_polygons);
	for (uint32_t i = 0; i < voronoi->N_polygons; i++)
		free(voronoi->vertices_forming_polygons[i]);
	free(voronoi->vertices_forming_polygons);
	if (voronoi->adjacencies != NULL) {
		for (uint32_t i = 0; i < voronoi->N_polygons; i++)
			free(voronoi->adjacencies[i]);
		free(voronoi->N_adjacencies);
		free(voronoi->adjacencies);
	}
	free(voronoi);
}

static inline int compare_voronoi_sgm(const void *const sgmA, 
				      const void *const sgmB)
{
	msh_edge_t* sA = ((voronoi_boundary_sgm_t*)sgmA)->sgm;
	msh_edge_t* sB = ((voronoi_boundary_sgm_t*)sgmB)->sgm;
	if (sA == sB)
		return 0;
	if (sA > sB)
		return 1;
	return -1;
}

static vcn_point2D_t* voronoi_insert_global(vcn_bins2D_t* bins, vcn_point2D_t* point)
{
	if(vcn_bins2D_are_points_inside_circle(bins, point->x, 
					       2.0 * NB_GEOMETRIC_TOL)) {
		vcn_point2D_t* point_out = NULL;
		double dist;
		vcn_bins2D_get_knn(bins, point, 1, &point_out, &dist);
		return point_out;
	}
	vcn_bins2D_insert(bins, point);
	return point;
}

static void voronoi_get_patched_circumcenter(double* v1,
					     double* v2,
					     double* v3,
					     /* Output */
					     double* circumcenter){
  /* Patch for bad quality triangles (Almost in boundaries) */
  double l1 = vcn_utils2D_get_dist(v1, v2);
  double l2 = vcn_utils2D_get_dist(v2, v3);
  double l3 = vcn_utils2D_get_dist(v3, v1);
  double* sv1 = v3;
  double* sv2 = v1;
  double cond;
  if(l1 > l2){
    if(l1 > l3){
      cond = l1 / (l2 + l3);
      sv1 = v1;
      sv2 = v2;
    }else
      cond = l3 / (l1 + l2);
  }else{
    if(l2 > l3){
      cond = l2 / (l1 + l3);
      sv1 = v2;
      sv2 = v3;
    }else
      cond = l3 / (l1 + l2);	      
  }
  if(cond > 0.8){
    /* Estimate the intersection made by the real
       circumcenter */
    circumcenter[0] = (sv1[0] + sv2[0]) * 0.5;
    circumcenter[1] = (sv1[1] + sv2[1]) * 0.5;	      
  }else{
    vcn_utils2D_get_circumcenter(v1, v2, v3,
		     circumcenter);
  }
}

static nb_container_t* voronoi_get_boundary_segments_with_vertices
(const vcn_mesh_t* const mesh)
{
	/* The size must not be zero if there are not input segments */
	nb_container_t* avl_boundary_sgm =
		nb_container_create(NB_SORTED);
	/* REFACTOR
	nb_container_set_key_generator(compare_voronoi_sgm);
	nb_container_set_comparer();
	*/
	for (uint32_t i=0; i < mesh->N_input_sgm; i++) {
		if (NULL == mesh->input_sgm[i])
			continue;
		if (NULL != mesh->input_sgm[i]->t1 && 
		    NULL != mesh->input_sgm[i]->t2)
			/* The input segment is not a boundary */
			continue;
		/* Compute boundary vertices of the segment */
		msh_edge_t* sgm = mesh->input_sgm[i];
		while (NULL != sgm) {
			voronoi_boundary_sgm_t* bsgm = calloc(1, sizeof(*bsgm));
			bsgm->sgm = sgm;
			bsgm->v1 = vcn_point2D_create();
			nb_container_insert(avl_boundary_sgm, bsgm);
			msh_trg_t* trg;
			if (sgm->t1 == NULL) {
				memcpy(bsgm->v1->x, sgm->v2->x,
				       2 * sizeof(double));
				trg = sgm->t2;
			} else {
				memcpy(bsgm->v1->x, sgm->v1->x,
				       2 * sizeof(double));
				trg = sgm->t1;
			}
			vcn_point2D_t* circumcenter = vcn_point2D_create();
			vcn_utils2D_get_circumcenter(trg->v1->x,
					 trg->v2->x,
					 trg->v3->x,
					 circumcenter->x);
			double area_test;
			if (sgm->t1 == trg)
				area_test = vcn_utils2D_get_2x_trg_area(sgm->v1->x,
							       sgm->v2->x,
							       circumcenter->x);
			else
				area_test = vcn_utils2D_get_2x_trg_area(sgm->v2->x,
							       sgm->v1->x,
							       circumcenter->x);
            
			if (area_test > NB_GEOMETRIC_TOL) {
				/* Case 1: The circumcenter is inside the boundary */
				bsgm->N_vertices = 1;
				circumcenter->x[0] = (sgm->v1->x[0] +
						      sgm->v2->x[0]) * 0.5;
				circumcenter->x[1] = (sgm->v1->x[1] +
						      sgm->v2->x[1]) * 0.5;
				bsgm->vertices[0] = circumcenter;
			} else if (area_test < 0) {
				/* Case 2: The circumcenter is outside the boundary */
				bsgm->N_vertices = 2;
				vcn_point2D_t* v1 = vcn_point2D_create();
				vcn_point2D_t* v2 = vcn_point2D_create();
				if (NULL == sgm->t2) {
					msh_vtx_t* sgm_v3 = 
						mtrg_get_opposite_vertex_guided
						(sgm->t1, sgm, true);
					double half1[2];
					half1[0] = (sgm->v1->x[0] +
						    sgm_v3->x[0])/2.0;
					half1[1] = (sgm->v1->x[1] +
						    sgm_v3->x[1])/2.0;
					vcn_utils2D_are_sgm_intersected
						(half1, circumcenter->x,
						 sgm->v1->x, sgm->v2->x,
						 v1->x, NULL);
					double half2[2];
					half2[0] = (sgm->v2->x[0] + 
						    sgm_v3->x[0])/2.0;
					half2[1] = (sgm->v2->x[1] +
						    sgm_v3->x[1])/2.0;
					vcn_utils2D_are_sgm_intersected
						(half2, circumcenter->x,
						 sgm->v1->x, sgm->v2->x,
						 v2->x, NULL);
				} else {
					msh_vtx_t* sgm_v3 = 
						mtrg_get_opposite_vertex_guided
						(sgm->t2, sgm, false);
					double half1[2];
					half1[0] = (sgm->v2->x[0] + 
						    sgm_v3->x[0])/2.0;
					half1[1] = (sgm->v2->x[1] +
						    sgm_v3->x[1])/2.0;
					vcn_utils2D_are_sgm_intersected
						(half1, circumcenter->x,
						 sgm->v1->x, sgm->v2->x,
						 v1->x, NULL);
					double half2[2];
					half2[0] = (sgm->v1->x[0] + 
						    sgm_v3->x[0])/2.0;
					half2[1] = (sgm->v1->x[1] +
						    sgm_v3->x[1])/2.0;
					vcn_utils2D_are_sgm_intersected
						(half2, circumcenter->x,
						 sgm->v1->x, sgm->v2->x,
						 v2->x, NULL);
				}
				bsgm->vertices[0] = v1;
				bsgm->vertices[1] = v2;
				free(circumcenter);
			} else {
				bsgm->N_vertices = 0;
			}
			sgm = medge_subsgm_next(sgm);
		}
	}
	return avl_boundary_sgm;
}


static int8_t voronoi_vtx_compare(const void *const A,
				  const void *const B)
{
	const double *const a = A;
	const double *const b = B;
	return ((fabs(a[0] - b[0]) < NB_GEOMETRIC_TOL &&
		 fabs(a[1] - b[1]) < NB_GEOMETRIC_TOL)) ? 1:0;
}

static void mesh_Lloyd_iteration(msh_vtx_t** vertices, uint32_t N_vertices,
				 msh_edge_t** segments, uint32_t N_segments,
				 msh_trg_t** triangles, uint32_t N_triangles,
				 /* NULL for a constant density */
				 double (*density)(double*),
				 const nb_container_t *const avl_boundary_sgm,
				 uint32_t max_iter)
{
	/* Initialize lists to store triangles connections */
	for (uint32_t i = 0; i < N_vertices; i++) {
		msh_vtx_t* vtx = vertices[i];
		void** attr = malloc(2 * sizeof(*attr));
		attr[0] = nb_container_create(NB_QUEUE);
		attr[1] = vtx->attr;
		vtx->attr = attr;
	}
	/* Get triangles connections per vertex */
	for (uint32_t i = 0; i < N_triangles; i++) {
		msh_trg_t* trg = triangles[i];
		/* Allocate circumcenter as an attribute */
		void** attr = malloc(2 * sizeof(*attr));
		double* circumcenter = malloc(2 * sizeof(*circumcenter));
		attr[0] = circumcenter;
		attr[1] = trg->attr;
		trg->attr = attr;
		/* Compute circumcenter */    
		voronoi_get_patched_circumcenter(trg->v1->x,
						 trg->v2->x,
						 trg->v3->x,
						 circumcenter);

		/* Set triangle to its vertices */
		nb_container_insert((nb_container_t*)((void**)trg->v1->attr)[0], trg);
		nb_container_insert((nb_container_t*)((void**)trg->v2->attr)[0], trg);
		nb_container_insert((nb_container_t*)((void**)trg->v3->attr)[0], trg);
	}
	/* Lloyd's iteration */
	double norm2 = 1e10;
	uint32_t k=0;
	while (true) {
		/* Move vertices to Voronoi centers of mass */
		norm2 = 0;
		for (uint32_t i=0; i < N_vertices; i++) {
			msh_vtx_t* vtx = vertices[i];
      
			/* Get vertices to calculate center of mass */
			nb_container_t* l_vtx = nb_container_create(NB_QUEUE);
			nb_container_set_comparer(l_vtx, voronoi_vtx_compare);
			nb_container_insert(l_vtx, vtx->x);
			nb_container_t* list = (nb_container_t*)((void**)vtx->attr)[0];
			nb_iterator_t* liter = nb_iterator_create();
			nb_iterator_set_container(liter, list);
			while (nb_iterator_has_more(liter)) {
				msh_trg_t* trg = (msh_trg_t*)nb_iterator_get_next(liter);
				double* circumcenter = (double*)((void**)trg->attr)[0];
	
				msh_edge_t* sgm_right = mtrg_get_right_edge(trg, vtx);
				msh_edge_t* sgm_left = mtrg_get_left_edge(trg, vtx);
				msh_edge_t* sgm_opposite = mtrg_get_opposite_edge(trg, vtx);
				if (medge_is_subsgm(sgm_opposite) &&
				    (sgm_opposite->t1 == NULL || sgm_opposite->t2 == NULL)) {
					voronoi_boundary_sgm_t aux_bsgm;
					aux_bsgm.sgm = sgm_opposite;
					voronoi_boundary_sgm_t* bsgm =
						nb_container_exist(avl_boundary_sgm, &aux_bsgm);
					if (bsgm->N_vertices == 2) {
						/* Case 2, opposite side */
						if (!nb_container_exist(l_vtx, bsgm->vertices[0]->x))
							nb_container_insert(l_vtx, bsgm->vertices[0]->x);
						if (!nb_container_exist(l_vtx, bsgm->vertices[1]->x))
							nb_container_insert(l_vtx, bsgm->vertices[1]->x);
					}
				}
	
				if (medge_is_subsgm(sgm_right) &&
				   (sgm_right->t1 == NULL || sgm_right->t2 == NULL)) {
					voronoi_boundary_sgm_t aux_bsgm;
					aux_bsgm.sgm = sgm_right;
					voronoi_boundary_sgm_t* bsgm =
						nb_container_exist(avl_boundary_sgm, &aux_bsgm);
					/* Add centroid to case 1 and 2 */
					if (!nb_container_exist(l_vtx, bsgm->v1->x))
						nb_container_insert(l_vtx, bsgm->v1->x);
					if (bsgm->N_vertices > 0) {
						if (!nb_container_exist(l_vtx, bsgm->vertices[0]->x))
							nb_container_insert(l_vtx, bsgm->vertices[0]->x);
					}
				}
				if (medge_is_subsgm(sgm_left) &&
				    (sgm_left->t1 == NULL || sgm_left->t2 == NULL)) {
					voronoi_boundary_sgm_t aux_bsgm;
					aux_bsgm.sgm = sgm_left;
					voronoi_boundary_sgm_t* bsgm =
						nb_container_exist(avl_boundary_sgm, &aux_bsgm);
	  
					/* Case 1 and 2 */
					if (bsgm->N_vertices > 0) {
						msh_vtx_t* bvtx = bsgm->vertices[0];
						if (bsgm->N_vertices == 2)
							bvtx = bsgm->vertices[1];
						if (!nb_container_exist(l_vtx, bvtx->x))
							nb_container_insert(l_vtx, bvtx->x);
					}
				}
	
				/* Always consider the circumcenter to compute the center of mass */
				if (!nb_container_exist(l_vtx, circumcenter))
					nb_container_insert(l_vtx, circumcenter);
			}
			nb_iterator_destroy(liter);

			/* Compute center of mass */
			double* internal_vtx =
				malloc(nb_container_get_length(l_vtx) * 2 * sizeof(*internal_vtx));
			uint32_t id = 0;
			liter = nb_iterator_create();
			nb_iterator_set_container(liter, l_vtx);
			while (nb_iterator_has_more(liter)) {
				double* p = (double*)nb_iterator_get_next(liter);
				memcpy(&(internal_vtx[id*2]), p, 2*sizeof(double));
				id++;
			}
			nb_iterator_destroy(liter);
			nb_container_destroy(l_vtx);
      			
			vcn_mesh_t* mesh_center_of_mass = vcn_mesh_create(); 
			vcn_mesh_get_delaunay(mesh_center_of_mass, id,
					      internal_vtx);

			free(internal_vtx);

			double next_x[2] = {0.0, 0.0};
			double normalizer = 0.0;
			nb_iterator_t* trg_iter = nb_iterator_create();
			nb_iterator_set_container(trg_iter, mesh_center_of_mass->ht_trg);
			while (nb_iterator_has_more(trg_iter)) {
				msh_trg_t* trg = (msh_trg_t*)nb_iterator_get_next(trg_iter);
				double trg_center[2];
				double trg_v1[2];
				trg_v1[0] = trg->v1->x[0]/mesh_center_of_mass->scale + mesh_center_of_mass->xdisp;
				trg_v1[1] = trg->v1->x[1]/mesh_center_of_mass->scale + mesh_center_of_mass->ydisp;
				double trg_v2[2];
				trg_v2[0] = trg->v2->x[0]/mesh_center_of_mass->scale + mesh_center_of_mass->xdisp;
				trg_v2[1] = trg->v2->x[1]/mesh_center_of_mass->scale + mesh_center_of_mass->ydisp;
				double trg_v3[2];
				trg_v3[0] = trg->v3->x[0]/mesh_center_of_mass->scale + mesh_center_of_mass->xdisp;
				trg_v3[1] = trg->v3->x[1]/mesh_center_of_mass->scale + mesh_center_of_mass->ydisp;
				if (density == NULL) {
					trg_center[0] = (trg_v1[0] + trg_v2[0] + trg_v3[0]) / 3.0;
					trg_center[1] = (trg_v1[1] + trg_v2[1] + trg_v3[1]) / 3.0;
				} else {
					double d1 = density(trg_v1);
					double d2 = density(trg_v2);
					double d3 = density(trg_v3);
					double density_sum = d1 + d2 + d3;
					trg_center[0] = 
						(trg_v1[0] * d1 + trg_v2[0] * d2 + trg_v3[0] * d3) / density_sum;
					trg_center[1] = 
						(trg_v1[1] * d1 + trg_v2[1] * d2 + trg_v3[1] * d3) / density_sum;
				}
	
				double area = 0.5 * vcn_utils2D_get_2x_trg_area(trg_v1, trg_v2, trg_v3);
				if (NULL == density) {
					next_x[0] += area * trg_center[0];
					next_x[1] += area * trg_center[1];
					normalizer += area;
				} else {
					double density_center = density(trg_center);
					next_x[0] += area * trg_center[0] * density_center;
					next_x[1] += area * trg_center[1] * density_center;
					normalizer += area * density_center;
				}
			}
			nb_iterator_destroy(trg_iter);
			vcn_mesh_destroy(mesh_center_of_mass);
      			
			next_x[0] /= normalizer;
			next_x[1] /= normalizer;

			/* Calculate norm for the stop criterion */
			norm2 += POW2(next_x[0] - vtx->x[0]) + 
				POW2(next_x[1] - vtx->x[1]);
      			
			/* Update position */
			memcpy(vtx->x, next_x, 2 * sizeof(*next_x));
		}
    
		/* Check stop criterion */
		k++;
		printf("min(voronoi): %e     (%i/%i) \r", norm2, k, max_iter); /* TEMPORAL */
		fflush(stdout);                                                /* TEMPORAL */
		if (norm2 < NB_GEOMETRIC_TOL)
			break;
		if (k > max_iter)
			break;

		/* Recompute circumcenters of triangles */
		for (uint32_t i = 0; i < N_triangles; i++) {
			msh_trg_t* trg = triangles[i];
			double center[2];
			vcn_utils2D_get_circumcenter(trg->v1->x,
					 trg->v2->x,
					 trg->v3->x,
					 center);
			double radii = vcn_utils2D_get_dist(trg->v1->x, center);
			bool flipped =  false;
			if (trg->t1 != NULL) {
				if (!medge_is_subsgm(trg->s1)) {
					msh_vtx_t* vtx = mtrg_get_opposite_vertex(trg->t1, trg->s1);
					if (radii - vcn_utils2D_get_dist(vtx->x, center) > NB_GEOMETRIC_TOL) {
						msh_edge_t* shared_sgm = trg->s1;
						msh_trg_t* t1 = shared_sgm->t1;
						msh_trg_t* t2 = shared_sgm->t2;
						nb_container_t* list1 = (nb_container_t*)((void**)shared_sgm->v1->attr)[0];
						nb_container_t* list2 = (nb_container_t*)((void**)shared_sgm->v2->attr)[0];
						nb_container_t* list3 = (nb_container_t*)
							((void**)mtrg_get_opposite_vertex(t1, shared_sgm)->attr)[0];
						nb_container_t* list4 = (nb_container_t*)
							((void**)mtrg_get_opposite_vertex(t2, shared_sgm)->attr)[0];
						/* Flip segment */
						medge_flip_without_dealloc(shared_sgm);
						/* Update triangular lists */
						nb_container_delete(list1, t1);
						nb_container_delete(list2, t2);
						nb_container_insert(list3, t2);
						nb_container_insert(list4, t1);

						flipped = true;
					}
				}
			}
			if (NULL != trg->t2 && !flipped){
				if (!medge_is_subsgm(trg->s2)) {
					msh_vtx_t* vtx = mtrg_get_opposite_vertex(trg->t2, trg->s2);
					if (radii - vcn_utils2D_get_dist(vtx->x, center) > NB_GEOMETRIC_TOL) {
						msh_edge_t* shared_sgm = trg->s2;
						msh_trg_t* t1 = shared_sgm->t1;
						msh_trg_t* t2 = shared_sgm->t2;
						nb_container_t* list1 = (nb_container_t*)((void**)shared_sgm->v1->attr)[0];
						nb_container_t* list2 = (nb_container_t*)((void**)shared_sgm->v2->attr)[0];
						nb_container_t* list3 = (nb_container_t*)
							((void**)mtrg_get_opposite_vertex(t1, shared_sgm)->attr)[0];
						nb_container_t* list4 = (nb_container_t*)
							((void**)mtrg_get_opposite_vertex(t2, shared_sgm)->attr)[0];
						/* Flip segment */
						medge_flip_without_dealloc(shared_sgm);
						/* Update triangular lists */
						nb_container_delete(list1, t1);
						nb_container_delete(list2, t2);
						nb_container_insert(list3, t2);
						nb_container_insert(list4, t1);

						flipped = true;
					}
				}
			}
			if (NULL != trg->t3 && !flipped) {
				if (!medge_is_subsgm(trg->s3)) {
					msh_vtx_t* vtx = mtrg_get_opposite_vertex(trg->t3, trg->s3);
					if (radii - vcn_utils2D_get_dist(vtx->x, center) > NB_GEOMETRIC_TOL) {
						msh_edge_t* shared_sgm = trg->s3;
						msh_trg_t* t1 = shared_sgm->t1;
						msh_trg_t* t2 = shared_sgm->t2;
						nb_container_t* list1 = (nb_container_t*)((void**)shared_sgm->v1->attr)[0];
						nb_container_t* list2 = (nb_container_t*)((void**)shared_sgm->v2->attr)[0];
						nb_container_t* list3 = (nb_container_t*)
							((void**)mtrg_get_opposite_vertex(t1, shared_sgm)->attr)[0];
						nb_container_t* list4 = (nb_container_t*)
							((void**)mtrg_get_opposite_vertex(t2, shared_sgm)->attr)[0];
						/* Flip segment */
						medge_flip_without_dealloc(shared_sgm);
						/* Update triangular lists */
						nb_container_delete(list1, t1);
						nb_container_delete(list2, t2);
						nb_container_insert(list3, t2);
						nb_container_insert(list4, t1);
					}
				}
			}
      
			double* circumcenter = (double*)((void**)trg->attr)[0];
			voronoi_get_patched_circumcenter(trg->v1->x,
							 trg->v2->x,
							 trg->v3->x,
							 circumcenter);
		}
	}
	printf("                                     \r"); fflush(stdout);/* Erase status */

	/* Free memory */
	for (uint32_t i = 0; i < N_vertices; i++) {
		msh_vtx_t* vtx = vertices[i];
		void** attr = (void**)vtx->attr;
		vtx->attr = attr[1];
		nb_container_destroy((nb_container_t*)attr[0]);
		free(attr);
	}

	for (uint32_t i = 0; i < N_triangles; i++) {
		msh_trg_t* trg = triangles[i];
		void** attr = (void**)trg->attr;
		trg->attr = attr[1];
		free(attr[0]);
		free(attr);
	}
}
