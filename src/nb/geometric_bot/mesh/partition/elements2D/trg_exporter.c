#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <alloca.h>

#include "nb/container_bot/container.h"
#include "nb/container_bot/iterator.h"
#include "nb/geometric_bot/knn/bins2D.h"
#include "nb/geometric_bot/knn/bins2D_iterator.h"
#include "nb/geometric_bot/mesh/partition/elements2D/trg_exporter.h"

#include "../../mesh2D_structs.h"

static void export_vertices(const vcn_mesh_t *const mesh,
			    nb_trg_exporter_interface_t *exp);
static void export_edges(const vcn_mesh_t *const mesh,
			 nb_trg_exporter_interface_t *exp);
static void export_and_enumerate_trg(vcn_mesh_t *mesh,
				     nb_trg_exporter_interface_t *exp);
static void export_trg_neighbours(const vcn_mesh_t *const mesh,
				  nb_trg_exporter_interface_t *exp);
static void export_input_vtx(const vcn_mesh_t *const mesh,
			     nb_trg_exporter_interface_t *exp);
static void export_input_sgm(const vcn_mesh_t *const mesh,
			     nb_trg_exporter_interface_t * exp);
static void set_input_sgm_table(const vcn_mesh_t *const restrict mesh,
				nb_trg_exporter_interface_t * restrict exp);
static void set_input_sgm(const vcn_mesh_t *const restrict mesh,
			  nb_trg_exporter_interface_t * restrict exp,
			  uint32_t isgm);

void vcn_mesh_export(const vcn_mesh_t *const mesh,
		     nb_trg_exporter_interface_t *exp)
{
	if (vcn_bins2D_is_empty(mesh->ug_vtx))
		goto EXIT;

	mesh_enumerate_vtx((vcn_mesh_t*)mesh);
	export_vertices(mesh, exp);

	if (NULL != exp->set_N_edg) {
		if (nb_container_is_not_empty(mesh->ht_edge))
			export_edges(mesh, exp);
	}

	if (NULL != exp->set_N_trg) {
		if (nb_container_is_not_empty(mesh->ht_trg)) {
			export_and_enumerate_trg((vcn_mesh_t*)mesh, exp);
			if (NULL != exp->set_trg_neighbours)
				export_trg_neighbours(mesh, exp);
		}
	}

	if (NULL != exp->set_N_input_vtx)
		export_input_vtx(mesh, exp);

	if (NULL != exp->set_N_input_sgm) {
		if (mesh->N_input_sgm > 0)
			export_input_sgm(mesh, exp);
	}

EXIT:
	return;
}

static void export_vertices(const vcn_mesh_t *const restrict mesh,
			    nb_trg_exporter_interface_t * restrict exp)
{
	uint32_t N_vtx = vcn_bins2D_get_length(mesh->ug_vtx);
	exp->set_N_vtx(exp->structure, N_vtx);
	exp->malloc_vtx(exp->structure);

	exp->start_vtx_access(exp->structure);

	vcn_bins2D_iter_t* iter = alloca(vcn_bins2D_iter_get_memsize());
	vcn_bins2D_iter_init(iter);
	vcn_bins2D_iter_set_bins(iter, mesh->ug_vtx);
	while (vcn_bins2D_iter_has_more(iter)) {
		const msh_vtx_t* vtx = vcn_bins2D_iter_get_next(iter);
		uint32_t id = mvtx_get_id(vtx);
		double x = vtx->x[0]/mesh->scale + mesh->xdisp;
		double y = vtx->x[1]/mesh->scale + mesh->ydisp;
		exp->set_vtx(exp->structure, id, x, y);
	}
	vcn_bins2D_iter_finish(iter);

	exp->stop_vtx_access(exp->structure);
}

static void export_edges(const vcn_mesh_t *const restrict mesh,
			 nb_trg_exporter_interface_t * restrict exp)
{
	uint32_t N_edg = nb_container_get_length(mesh->ht_edge);
	exp->set_N_edg(exp->structure, N_edg);
	exp->malloc_edg(exp->structure);

	exp->start_edg_access(exp->structure);

	uint32_t i = 0;
	nb_iterator_t *iter = alloca(nb_iterator_get_memsize());
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, mesh->ht_edge);
	while (nb_iterator_has_more(iter)) {
		msh_edge_t *edge = (msh_edge_t*) nb_iterator_get_next(iter);
		uint32_t v1 = mvtx_get_id(edge->v1);
		uint32_t v2 = mvtx_get_id(edge->v2);
		exp->set_edg(exp->structure, i, v1, v2);
		i += 1;
	}
	nb_iterator_finish(iter);

	exp->stop_edg_access(exp->structure);
}

static void export_and_enumerate_trg(vcn_mesh_t *mesh,
				    nb_trg_exporter_interface_t *exp)
{
	uint32_t N_trg = nb_container_get_length(mesh->ht_trg);
	exp->set_N_trg(exp->structure, N_trg);
	bool include_neighbours = (NULL != exp->set_trg_neighbours);
	exp->malloc_trg(exp->structure, include_neighbours);

	exp->start_trg_access(exp->structure);

	nb_iterator_t* trg_iter = alloca(nb_iterator_get_memsize());
	nb_iterator_init(trg_iter);
	nb_iterator_set_container(trg_iter, mesh->ht_trg);
	uint32_t i = 0;
	while (nb_iterator_has_more(trg_iter)) {
		msh_trg_t* trg = (msh_trg_t*)nb_iterator_get_next(trg_iter);
		uint32_t v1 = mvtx_get_id(trg->v1);
		uint32_t v2 = mvtx_get_id(trg->v2);
		uint32_t v3 = mvtx_get_id(trg->v3);
		exp->set_trg(exp->structure, i, v1, v2, v3);

		trg->id = i;
		i++;
	}
	nb_iterator_finish(trg_iter);

	exp->stop_trg_access(exp->structure);
}

static void export_trg_neighbours(const vcn_mesh_t *const restrict mesh,
				  nb_trg_exporter_interface_t * restrict exp)
{
	exp->start_trg_neighbours_access(exp->structure);

	nb_iterator_t* trg_iter = alloca(nb_iterator_get_memsize());
	nb_iterator_init(trg_iter);
	nb_iterator_set_container(trg_iter, mesh->ht_trg);
	while (nb_iterator_has_more(trg_iter)) {
		msh_trg_t* trg = (msh_trg_t*)nb_iterator_get_next(trg_iter);
		uint32_t id = trg->id;
		uint32_t t1 = vcn_mesh_get_N_trg(mesh);
		if (NULL != trg->t1)
			t1 = trg->t1->id;
		uint32_t t2 = vcn_mesh_get_N_trg(mesh);
		if (NULL != trg->t2)
			t2 = trg->t2->id;
		uint32_t t3 = vcn_mesh_get_N_trg(mesh);
		if (NULL != trg->t3)
			t3 = trg->t3->id;

		exp->set_trg_neighbours(exp->structure, id, t1, t2, t3);
	}
	nb_iterator_finish(trg_iter);

	exp->stop_trg_neighbours_access(exp->structure);
}

static void export_input_vtx(const vcn_mesh_t *const restrict mesh,
			     nb_trg_exporter_interface_t * restrict exp)
{
	uint32_t N_vtx = mesh->N_input_vtx;
	exp->set_N_input_vtx(exp->structure, N_vtx);
	exp->malloc_input_vtx(exp->structure);

	exp->start_input_vtx_access(exp->structure);
	for (uint32_t i = 0; i < N_vtx; i++) {
		uint32_t vtx_id;
		if (NULL == mesh->input_vtx[i])
			vtx_id = vcn_mesh_get_N_vtx(mesh);
		else
			vtx_id = mvtx_get_id(mesh->input_vtx[i]);
		exp->set_input_vtx(exp->structure, i, vtx_id);
	}
	exp->stop_input_vtx_access(exp->structure);
}

static void export_input_sgm(const vcn_mesh_t *const restrict mesh,
			     nb_trg_exporter_interface_t * restrict exp)
{
	uint32_t N_sgm = mesh->N_input_sgm;
	exp->set_N_input_sgm(exp->structure, N_sgm);
	exp->malloc_input_sgm_table(exp->structure);

	exp->start_input_sgm_table_access(exp->structure);
	set_input_sgm_table(mesh, exp);
	
	for (uint32_t i = 0; i < N_sgm; i++) {
		if (NULL != mesh->input_sgm[i]) {
			exp->input_sgm_start_access(exp->structure, i);
			set_input_sgm(mesh, exp, i);
			exp->input_sgm_stop_access(exp->structure);
		}
	}
	exp->stop_input_sgm_table_access(exp->structure);
}


static void set_input_sgm_table(const vcn_mesh_t *const restrict mesh,
				nb_trg_exporter_interface_t * restrict exp)
{
	for (uint32_t i = 0; i < mesh->N_input_sgm; i++) {
		msh_edge_t* sgm = mesh->input_sgm[i];
		uint32_t counter = 0;
		while (NULL != sgm) {
			counter++;
			sgm = medge_subsgm_next(sgm);
		}
		if (0 < counter) {
			exp->input_sgm_set_N_vtx(exp->structure, i, counter + 1);
			exp->input_sgm_malloc_vtx(exp->structure, i);
		} else {
			exp->input_sgm_set_N_vtx(exp->structure, i, 0);
		}
	}
}

static void set_input_sgm(const vcn_mesh_t *const restrict mesh,
			  nb_trg_exporter_interface_t * restrict exp,
			  uint32_t isgm)
{
	msh_edge_t* sgm = mesh->input_sgm[isgm];
	msh_edge_t* sgm_prev = sgm;
	sgm = medge_subsgm_next(sgm);
	if (NULL == sgm) {
		exp->input_sgm_set_vtx(exp->structure, 0,
				       mvtx_get_id(sgm_prev->v1));
		exp->input_sgm_set_vtx(exp->structure, 1,
				       mvtx_get_id(sgm_prev->v2));
	} else {
		uint32_t idx = 0;
		uint32_t id_chain;
		uint32_t id1 = mvtx_get_id(sgm_prev->v1);
		uint32_t id2 = mvtx_get_id(sgm_prev->v2);
		uint32_t id1n = mvtx_get_id(sgm->v1);
		uint32_t id2n = mvtx_get_id(sgm->v2);
		if (id2 == id1n || id2 == id2n) {
			exp->input_sgm_set_vtx(exp->structure, idx++, id1);
			exp->input_sgm_set_vtx(exp->structure, idx++, id2);
			id_chain = id2;
		} else {
			exp->input_sgm_set_vtx(exp->structure, idx++, id2);
			exp->input_sgm_set_vtx(exp->structure, idx++, id1);
			id_chain = id1;
		}
		while (NULL != sgm) {
			sgm_prev = sgm;
			uint32_t id1 = mvtx_get_id(sgm_prev->v1);
			uint32_t id2 = mvtx_get_id(sgm_prev->v2);
			if (id1 == id_chain) {
				exp->input_sgm_set_vtx(exp->structure,
						       idx++, id2);
				id_chain = id2;
			} else {
				exp->input_sgm_set_vtx(exp->structure,
						       idx++, id1);
				id_chain = id1;
			}
			sgm = medge_subsgm_next(sgm);
		}
	}
}
