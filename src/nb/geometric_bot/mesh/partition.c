#include <stdlib.h>
#include <stdint.h>

#include "nb/graph_bot.h"
#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/geometric_bot/mesh/elements2D/triangles.h"
#include "nb/geometric_bot/mesh/elements2D/polygons.h"
#include "nb/geometric_bot/mesh/elements2D/quad.h"
#include "nb/geometric_bot/mesh/elements2D/disks.h"

#include "nb/geometric_bot/mesh/partition.h"

struct nb_partition_s {
	void *msh;
	nb_partition_type type;
	void (*finish)(void *msh);
	void (*copy)(void *msh, const void *mshsrc);
	void (*clear)(void *msh);
	uint32_t (*get_N_invtx)(const void *msh);
	uint32_t (*get_N_insgm)(const void *msh);
	uint32_t (*get_N_nodes)(const void *msh);
	uint32_t (*get_N_edges)(const void *msh);
	uint32_t (*get_N_elems)(const void *msh);
	double (*get_x_node)(const void *msh, uint32_t id);
	double (*get_y_node)(const void *msh, uint32_t id);
	uint32_t (*get_1n_edge)(const void *msh, uint32_t id);
	uint32_t (*get_2n_edge)(const void *msh, uint32_t id);
	double (*get_x_elem)(const void *msh, uint32_t id);
	double (*get_y_elem)(const void *msh, uint32_t id);
	uint32_t (*elem_get_N_adj)(const void *msh, uint32_t id);
	uint32_t (*elem_get_adj)(const void *msh,
				 uint32_t elem_id, uint8_t adj_id);
	uint32_t (*elem_get_N_ngb)(const void *msh, uint32_t id);
	uint32_t (*elem_get_ngb)(const void *msh,
				 uint32_t elem_id, uint8_t ngb_id);
	uint32_t (*get_invtx)(const void *msh, uint32_t id);
	uint32_t (*get_N_nodes_x_insgm)(const void *msh, uint32_t id);
	uint32_t (*get_node_x_insgm)(const void *msh, uint32_t sgm_id,
				     uint32_t node_id);
	void (*load_elem_graph)(const void *msh, nb_graph_t *graph);
	void (*load_nodal_graph)(const void *msh, nb_graph_t *graph);
	void (*load_interelem_graph)(const void *msh, nb_graph_t *graph);
	void (*load_from_mesh)(void *msh, const nb_mesh_t *mesh);
	void (*get_enveloping_box)(const void *msh, double box[4]);
	bool (*is_vtx_inside)(const void *msh, double x, double y);
	void (*draw)(const void *msh, const char *filename,
		     int width, int height);
	void (*build_model)(const void *msh, nb_model_t *model);
	void (*build_model_disabled_elems)(const void *msh,
					   const bool *elems_enabled,
					   nb_model_t *model,
					   uint32_t *N_input_vtx,
					   uint32_t **input_vtx);
};

static void set_msh3trg_interface(nb_partition_t *part);
static void set_mshquad_interface(nb_partition_t *part);
static void set_mshpoly_interface(nb_partition_t *part);
static void set_mshpack_interface(nb_partition_t *part);

uint32_t nb_partition_get_memsize(nb_partition_type  type)
{
	uint32_t mem;
	switch (type) {
	case NB_TRIAN:
		mem = nb_msh3trg_get_memsize();
		break;
	case NB_QUAD:
		mem = nb_mshquad_get_memsize();
		break;
	case NB_POLY:
		mem = nb_mshpoly_get_memsize();
		break;
	case NB_DISK:
		mem = nb_mshpack_get_memsize();
		break;
	default:
		mem = nb_msh3trg_get_memsize();
		break;
	}
	return mem + sizeof(nb_partition_t);
}

void nb_partition_init(nb_partition_t *part, nb_partition_type  type)
{
	char *memblock = (void*) part;
	part->msh = (void*) (memblock + sizeof(nb_partition_t));
	part->type = type;
	switch (type) {
	case NB_TRIAN:
		nb_msh3trg_init(part->msh);
		set_msh3trg_interface(part);
		break;
	case NB_QUAD:
		nb_mshquad_init(part->msh);
		set_mshquad_interface(part);
		break;
	case NB_POLY:
		nb_mshpoly_init(part->msh);
		set_mshpoly_interface(part);
		break;
	case NB_DISK:
		nb_mshpack_init(part->msh);
		set_mshpack_interface(part);
		break;
	default:
		nb_msh3trg_init(part->msh);
		set_msh3trg_interface(part);
		break;
	}
}

static void set_msh3trg_interface(nb_partition_t *part)
{
	part->finish = nb_msh3trg_finish;
	part->copy = nb_msh3trg_copy;
	part->clear = nb_msh3trg_clear;
	part->get_N_invtx = nb_msh3trg_get_N_invtx;
	part->get_N_insgm = nb_msh3trg_get_N_insgm;
	part->get_N_nodes = nb_msh3trg_get_N_nodes;
	part->get_N_edges = nb_msh3trg_get_N_edges;
	part->get_N_elems = nb_msh3trg_get_N_elems;
	part->get_x_node = nb_msh3trg_get_x_node;
	part->get_y_node = nb_msh3trg_get_y_node;
	part->get_1n_edge = nb_msh3trg_get_1n_edge;
	part->get_2n_edge = nb_msh3trg_get_2n_edge;
	part->get_x_elem = nb_msh3trg_get_x_elem;
	part->get_y_elem = nb_msh3trg_get_y_elem;
	part->elem_get_N_adj = nb_msh3trg_elem_get_N_adj;
	part->elem_get_adj = nb_msh3trg_elem_get_adj;
	part->elem_get_N_ngb = nb_msh3trg_elem_get_N_ngb;
	part->elem_get_ngb = nb_msh3trg_elem_get_ngb;
	part->get_invtx = nb_msh3trg_get_invtx;
	part->get_N_nodes_x_insgm = nb_msh3trg_get_N_nodes_x_insgm;
	part->get_node_x_insgm = nb_msh3trg_get_node_x_insgm;
	part->load_elem_graph = nb_msh3trg_load_elem_graph;
	part->load_nodal_graph = nb_msh3trg_load_nodal_graph;
	part->load_interelem_graph = nb_msh3trg_load_interelem_graph;
	part->load_from_mesh = nb_msh3trg_load_from_mesh;
	part->get_enveloping_box = nb_msh3trg_get_enveloping_box;
	part->is_vtx_inside = nb_msh3trg_is_vtx_inside;
	part->draw = nb_msh3trg_draw;
	part->build_model = nb_msh3trg_build_model;
	part->build_model_disabled_elems =
		nb_msh3trg_build_model_disabled_elems;
}

static void set_mshquad_interface(nb_partition_t *part)
{
	part->finish = nb_mshquad_finish;
	part->copy = nb_mshquad_copy;
	part->clear = nb_mshquad_clear;
	part->get_N_invtx = nb_mshquad_get_N_invtx;
	part->get_N_insgm = nb_mshquad_get_N_insgm;
	part->get_N_nodes = nb_mshquad_get_N_nodes;
	part->get_N_edges = nb_mshquad_get_N_edges;
	part->get_N_elems = nb_mshquad_get_N_elems;
	part->get_x_node = nb_mshquad_get_x_node;
	part->get_y_node = nb_mshquad_get_y_node;
	part->get_1n_edge = nb_mshquad_get_1n_edge;
	part->get_2n_edge = nb_mshquad_get_2n_edge;
	part->get_x_elem = nb_mshquad_get_x_elem;
	part->get_y_elem = nb_mshquad_get_y_elem;
	part->elem_get_N_adj = nb_mshquad_elem_get_N_adj;
	part->elem_get_adj = nb_mshquad_elem_get_adj;
	part->elem_get_N_ngb = nb_mshquad_elem_get_N_ngb;
	part->elem_get_ngb = nb_mshquad_elem_get_ngb;
	part->get_invtx = nb_mshquad_get_invtx;
	part->get_N_nodes_x_insgm = nb_mshquad_get_N_nodes_x_insgm;
	part->get_node_x_insgm = nb_mshquad_get_node_x_insgm;
	part->load_elem_graph = nb_mshquad_load_elem_graph;
	part->load_nodal_graph = nb_mshquad_load_nodal_graph;
	part->load_interelem_graph = nb_mshquad_load_interelem_graph;
	part->load_from_mesh = nb_mshquad_load_from_mesh;
	part->get_enveloping_box = nb_mshquad_get_enveloping_box;
	part->is_vtx_inside = nb_mshquad_is_vtx_inside;
	part->draw = nb_mshquad_draw;
	part->build_model = nb_mshquad_build_model;
	part->build_model_disabled_elems =
		nb_mshquad_build_model_disabled_elems;
}

static void set_mshpoly_interface(nb_partition_t *part)
{
	part->finish = nb_mshpoly_finish;
	part->copy = nb_mshpoly_copy;
	part->clear = nb_mshpoly_clear;
	part->get_N_invtx = nb_mshpoly_get_N_invtx;
	part->get_N_insgm = nb_mshpoly_get_N_insgm;
	part->get_N_nodes = nb_mshpoly_get_N_nodes;
	part->get_N_edges = nb_mshpoly_get_N_edges;
	part->get_N_elems = nb_mshpoly_get_N_elems;
	part->get_x_node = nb_mshpoly_get_x_node;
	part->get_y_node = nb_mshpoly_get_y_node;
	part->get_1n_edge = nb_mshpoly_get_1n_edge;
	part->get_2n_edge = nb_mshpoly_get_2n_edge;
	part->get_x_elem = nb_mshpoly_get_x_elem;
	part->get_y_elem = nb_mshpoly_get_y_elem;
	part->elem_get_N_adj = nb_mshpoly_elem_get_N_adj;
	part->elem_get_adj = nb_mshpoly_elem_get_adj;
	part->elem_get_N_ngb = nb_mshpoly_elem_get_N_ngb;
	part->elem_get_ngb = nb_mshpoly_elem_get_ngb;
	part->get_invtx = nb_mshpoly_get_invtx;
	part->get_N_nodes_x_insgm = nb_mshpoly_get_N_nodes_x_insgm;
	part->get_node_x_insgm = nb_mshpoly_get_node_x_insgm;
	part->load_elem_graph = nb_mshpoly_load_elem_graph;
	part->load_nodal_graph = nb_mshpoly_load_nodal_graph;
	part->load_interelem_graph = nb_mshpoly_load_interelem_graph;
	part->load_from_mesh = nb_mshpoly_load_from_mesh;
	part->get_enveloping_box = nb_mshpoly_get_enveloping_box;
	part->is_vtx_inside = nb_mshpoly_is_vtx_inside;
	part->draw = nb_mshpoly_draw;
	part->build_model = nb_mshpoly_build_model;
	part->build_model_disabled_elems =
		nb_mshpoly_build_model_disabled_elems;
}

static void set_mshpack_interface(nb_partition_t *part)
{
	part->finish = nb_mshpack_finish;
	part->copy = nb_mshpack_copy;
	part->clear = nb_mshpack_clear;
	part->get_N_invtx = nb_mshpack_get_N_invtx;
	part->get_N_insgm = nb_mshpack_get_N_insgm;
	part->get_N_nodes = nb_mshpack_get_N_nodes;
	part->get_N_edges = nb_mshpack_get_N_edges;
	part->get_N_elems = nb_mshpack_get_N_elems;
	part->get_x_node = nb_mshpack_get_x_node;
	part->get_y_node = nb_mshpack_get_y_node;
	part->get_1n_edge = nb_mshpack_get_1n_edge;
	part->get_2n_edge = nb_mshpack_get_2n_edge;
	part->get_x_elem = nb_mshpack_get_x_elem;
	part->get_y_elem = nb_mshpack_get_y_elem;
	part->elem_get_N_adj = nb_mshpack_elem_get_N_adj;
	part->elem_get_adj = nb_mshpack_elem_get_adj;
	part->elem_get_N_ngb = nb_mshpack_elem_get_N_ngb;
	part->elem_get_ngb = nb_mshpack_elem_get_ngb;
	part->get_invtx = nb_mshpack_get_invtx;
	part->get_N_nodes_x_insgm = nb_mshpack_get_N_nodes_x_insgm;
	part->get_node_x_insgm = nb_mshpack_get_node_x_insgm;
	part->load_elem_graph = nb_mshpack_load_elem_graph;
	part->load_nodal_graph = nb_mshpack_load_nodal_graph;
	part->load_interelem_graph = nb_mshpack_load_interelem_graph;
	part->load_from_mesh = nb_mshpack_load_from_mesh;
	part->get_enveloping_box = nb_mshpack_get_enveloping_box;
	part->is_vtx_inside = nb_mshpack_is_vtx_inside;
	part->draw = nb_mshpack_draw;
	part->build_model = nb_mshpack_build_model;
	part->build_model_disabled_elems =
		nb_mshpack_build_model_disabled_elems;
}

void nb_partition_copy(nb_partition_t *part, const nb_partition_t* srcpart)
{
	nb_partition_init(part, srcpart->type);
	part->copy(part->msh, srcpart->msh);
}

void nb_partition_finish(nb_partition_t *part)
{
	part->finish(part->msh);
}

nb_partition_t* nb_partition_create(nb_partition_type type)
{
	uint32_t memsize = nb_partition_get_memsize(type);
	nb_partition_t *part = malloc(memsize);
	nb_partition_init(part, type);
	return part;
}

nb_partition_t* nb_partition_clone(nb_partition_t* part)
{
	uint32_t memsize = nb_partition_get_memsize(part->type);
	nb_partition_t *clone = malloc(memsize);
	nb_partition_copy(clone, part);
	return clone;
}

void nb_partition_clear(nb_partition_t* part)
{
	part->clear(part->msh);
}

void nb_partition_destroy(nb_partition_t* part)
{
	nb_partition_finish(part);
	free(part);
}

nb_partition_type nb_partition_get_type(const nb_partition_t *part)
{
	return part->type;
}

uint32_t nb_partition_get_N_invtx(const nb_partition_t *part)
{
	return part->get_N_invtx(part->msh);
}

uint32_t nb_partition_get_N_insgm(const nb_partition_t *part)
{
	return part->get_N_insgm(part->msh);
}

uint32_t nb_partition_get_N_nodes(const nb_partition_t *part)
{
	return part->get_N_nodes(part->msh);
}

uint32_t nb_partition_get_N_edges(const nb_partition_t *part)
{
	return part->get_N_edges(part->msh);
}

uint32_t nb_partition_get_N_elems(const nb_partition_t *part)
{
	return part->get_N_elems(part->msh);
}

double nb_partition_get_x_node(const nb_partition_t *part, uint32_t id)
{
	return part->get_x_node(part->msh, id);
}

double nb_partition_get_y_node(const nb_partition_t *part, uint32_t id)
{
	return part->get_y_node(part->msh, id);
}

uint32_t nb_partition_get_1n_edge(const nb_partition_t *part, uint32_t id)
{
	return part->get_1n_edge(part->msh, id);
}

uint32_t nb_partition_get_2n_edge(const nb_partition_t *part, uint32_t id)
{
	return part->get_2n_edge(part->msh, id);
}

double nb_partition_get_x_elem(const nb_partition_t *part, uint32_t id)
{
	return part->get_x_elem(part->msh, id);
}

double nb_partition_get_y_elem(const nb_partition_t *part, uint32_t id)
{
	return part->get_y_elem(part->msh, id);
}

uint32_t nb_partition_elem_get_N_adj(const nb_partition_t *part, uint32_t id)
{
	return part->elem_get_N_adj(part->msh, id);
}

uint32_t nb_partition_elem_get_adj(const nb_partition_t *part,
				   uint32_t elem_id, uint8_t adj_id)
{
	return part->elem_get_adj(part->msh, elem_id, adj_id);
}

uint32_t nb_partition_elem_get_N_ngb(const nb_partition_t *part, uint32_t id)
{
	return part->elem_get_N_ngb(part->msh, id);
}

uint32_t nb_partition_elem_get_ngb(const nb_partition_t *part,
				   uint32_t elem_id, uint8_t ngb_id)
{
	return part->elem_get_ngb(part->msh, elem_id, ngb_id);
}

uint32_t nb_partition_get_invtx(const nb_partition_t *part, uint32_t id)
{
	return part->get_invtx(part->msh, id);
}

uint32_t nb_partition_get_N_nodes_x_insgm(const nb_partition_t *part,
					  uint32_t id)
{
	return part->get_N_nodes_x_insgm(part->msh, id);
}

uint32_t nb_partition_get_node_x_insgm(const nb_partition_t *part,
				       uint32_t sgm_id, uint32_t node_id)
{
	return part->get_node_x_insgm(part->msh, sgm_id, node_id);
}

void nb_partition_load_elem_graph(const nb_partition_t *part,
				  nb_graph_t *graph)
{
	part->load_elem_graph(part->msh, graph);
}

void nb_partition_load_nodal_graph(const nb_partition_t *part,
				   vcn_graph_t *graph)
{
	part->load_nodal_graph(part->msh, graph);
}

void nb_partition_load_interelem_graph(const nb_partition_t *part,
				       vcn_graph_t *graph)
{
	part->load_interelem_graph(part->msh, graph);
}

void nb_partition_load_from_mesh(nb_partition_t *part,
				 const nb_mesh_t *const mesh)
{
	part->load_from_mesh(part->msh, mesh);
}

void nb_partition_get_enveloping_box(const nb_partition_t *part,
				     double box[4])
{
	part->get_enveloping_box(part->msh, box);
}

bool nb_partition_is_vtx_inside(const nb_partition_t *part,
				double x, double y)
{
	return part->is_vtx_inside(part->msh, x, y);
}

void nb_partition_draw(const nb_partition_t *part, const char *filename,
		       int width, int height)
{
	part->draw(part->msh, filename, width, height);
}

void nb_partition_build_model(const nb_partition_t *part, nb_model_t *model)
{
	part->build_model(part->msh, model);
}

void nb_partition_build_model_disabled_elems
			(const nb_partition_t *part,
			 const bool *elems_enabled,
			 nb_model_t *model,
			 uint32_t *N_input_vtx,
			 uint32_t **input_vtx)
{
	part->build_model_disabled_elems(part->msh,
					 elems_enabled,
					 model,
					 N_input_vtx,
					 input_vtx);
}
