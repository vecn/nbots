#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "nb/container_bot.h"
#include "nb/geometric_bot/knn/bins2D.h"
#include "nb/geometric_bot/knn/bins2D_iterator.h"
#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/geometric_bot/mesh/modules2D/extra_fem.h"

#include "../mesh2D_structs.h"

void vcn_mesh_duplicate_one_point_connections(vcn_mesh_t* mesh)
{
	/* Allocate lists to store triangles per vertex */
	vcn_bins2D_iter_t* iter = vcn_bins2D_iter_create();
	vcn_bins2D_iter_set_bins(iter, mesh->ug_vtx);
	while (vcn_bins2D_iter_has_more(iter)) {
		msh_vtx_t* vtx = (msh_vtx_t*) vcn_bins2D_iter_get_next(iter);
		void** attr = malloc(2 * sizeof(*attr));
		nb_container_t* trg_x_vtx =
			nb_container_create(NB_QUEUE);
		attr[0] = trg_x_vtx;
		attr[1] = vtx->attr;
		vtx->attr = attr;
	}
	/* Iterate over triangles to found relations */
	nb_iterator_t* ht_iter = nb_iterator_create();
	nb_iterator_set_container(ht_iter, mesh->ht_trg);
	while (nb_iterator_has_more(ht_iter)) {
		const msh_trg_t* trg = nb_iterator_get_next(ht_iter);
		nb_container_insert(((void**)trg->v1->attr)[0], trg);
		nb_container_insert(((void**)trg->v2->attr)[0], trg);
		nb_container_insert(((void**)trg->v3->attr)[0], trg);
	}
	nb_iterator_destroy(ht_iter);

	/* Detect one point connections and duplicate vertices */
	nb_container_t* new_vertices =
		nb_container_create(NB_QUEUE);
	vcn_bins2D_iter_restart(iter);
	while (vcn_bins2D_iter_has_more(iter)) {
		const msh_vtx_t* vtx = vcn_bins2D_iter_get_next(iter);
		nb_container_t* trg_x_vtx = ((void**)vtx->attr)[0];
		if (2 > nb_container_get_length(trg_x_vtx))
			continue;
		msh_trg_t* trg = nb_container_get_first(trg_x_vtx);

		msh_trg_t* trg_twist = trg;
		bool twist_around = false;
		while (NULL != mtrg_get_right_triangle(trg_twist, vtx)) {
			trg_twist = mtrg_get_right_triangle(trg_twist, vtx);
			if (trg_twist == trg) {
				twist_around = true;
				break;
			}
		}
		if (trg_twist == trg && twist_around)
			continue;

		nb_container_t* trg_fan =
			nb_container_create(NB_QUEUE);
		do {
			nb_container_insert(trg_fan, trg_twist);
			trg_twist = mtrg_get_left_triangle(trg_twist, vtx);
		} while(NULL != trg_twist);

		if (nb_container_get_length(trg_fan) == 
		    nb_container_get_length(trg_x_vtx)) {
			nb_container_destroy(trg_fan);
			continue;
		}
    
		msh_vtx_t* new_vtx = mvtx_create(mesh);
		memcpy(new_vtx->x, vtx->x, 2 * sizeof(*(vtx->x)));
		new_vtx->attr = ((void**)vtx->attr)[1];
		nb_container_insert(new_vertices, new_vtx);
		
		while (nb_container_is_not_empty(trg_fan)) {
			trg = nb_container_delete_first(trg_fan);
    
			msh_edge_t* s1 = mtrg_get_left_edge(trg, vtx);
			msh_edge_t* s2 = mtrg_get_right_edge(trg, vtx);
    
			if (trg->v1 == vtx)
				trg->v1 = new_vtx;
			else if (trg->v2 == vtx)
				trg->v2 = new_vtx;
			else if (trg->v3 == vtx)
				trg->v3 = new_vtx;
    
			if (s1->v1 == vtx)
				s1->v1 = new_vtx;
			else if (s1->v2 == vtx)
				s1->v2 = new_vtx;
    
			if (s2->v1 == vtx)
				s2->v1 = new_vtx;
			else if (s2->v2 == vtx)
				s2->v2 = new_vtx;
		}
    		nb_container_destroy(trg_fan);
	}

	/* Free memory */
	vcn_bins2D_iter_restart(iter);
	while (vcn_bins2D_iter_has_more(iter)) {
		msh_vtx_t* vtx = (msh_vtx_t*) vcn_bins2D_iter_get_next(iter);
		void** attr = vtx->attr;
		nb_container_destroy(attr[0]);
		vtx->attr = attr[1];
		free(attr);
	}
	vcn_bins2D_iter_destroy(iter);

	while (nb_container_is_not_empty(new_vertices)) {
		msh_vtx_t* new_vtx = nb_container_delete_first(new_vertices);
		vcn_bins2D_insert(mesh->ug_vtx, new_vtx);
	}

	nb_container_destroy(new_vertices);
}
