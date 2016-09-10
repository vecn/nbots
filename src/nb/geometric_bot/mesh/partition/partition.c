#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "nb/graph_bot.h"
#include "nb/graphics_bot.h"
#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/geometric_bot/mesh/partition/elements2D/msh3trg.h"
#include "nb/geometric_bot/mesh/partition/elements2D/msh3trg_draw.h"
#include "nb/geometric_bot/mesh/partition/elements2D/mshquad.h"
#include "nb/geometric_bot/mesh/partition/elements2D/mshquad_draw.h"
#include "nb/geometric_bot/mesh/partition/elements2D/mshpoly.h"
#include "nb/geometric_bot/mesh/partition/elements2D/mshpoly_draw.h"
#include "nb/geometric_bot/mesh/partition/elements2D/mshpack.h"
#include "nb/geometric_bot/mesh/partition/elements2D/mshpack_draw.h"

#include "nb/geometric_bot/mesh/partition.h"
#include "nb/geometric_bot/mesh/partition/info.h"
#include "partition_struct.h"

#define POW2(a) ((a)*(a))

static void set_msh_interface(nb_partition_t *part, nb_partition_type  type);

static void set_msh3trg_interface(nb_partition_t *part);
static void set_msh3trg_main_interface(nb_partition_t *part);
static void set_msh3trg_graphics_interface(nb_partition_t *part);

static void set_mshquad_interface(nb_partition_t *part);
static void set_mshquad_main_interface(nb_partition_t *part);
static void set_mshquad_graphics_interface(nb_partition_t *part);

static void set_mshpoly_interface(nb_partition_t *part);
static void set_mshpoly_main_interface(nb_partition_t *part);
static void set_mshpoly_graphics_interface(nb_partition_t *part);

static void set_mshpack_interface(nb_partition_t *part);
static void set_mshpack_main_interface(nb_partition_t *part);
static void set_mshpack_graphics_interface(nb_partition_t *part);

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
	set_msh_interface(part, type);
	part->init(part->msh);
}

static void set_msh_interface(nb_partition_t *part, nb_partition_type  type)
{
	switch (type) {
	case NB_TRIAN:
		set_msh3trg_interface(part);
		break;
	case NB_QUAD:
		set_mshquad_interface(part);
		break;
	case NB_POLY:
		set_mshpoly_interface(part);
		break;
	case NB_DISK:
		set_mshpack_interface(part);
		break;
	default:
		set_msh3trg_interface(part);
		break;
	}
}

void nb_partition_init_from_msh(nb_partition_t *part, void *msh,
				nb_partition_type  type)
{
	char *memblock = (void*) part;
	part->msh = msh;
	part->type = type;
	set_msh_interface(part, type);
}

static void set_msh3trg_interface(nb_partition_t *part)
{
	set_msh3trg_main_interface(part);
	set_msh3trg_graphics_interface(part);
}

static void set_msh3trg_main_interface(nb_partition_t *part)
{
	part->init = nb_msh3trg_init;
	part->finish = nb_msh3trg_finish;
	part->copy = nb_msh3trg_copy;
	part->clear = nb_msh3trg_clear;
	part->get_N_invtx = nb_msh3trg_get_N_invtx;
	part->get_N_insgm = nb_msh3trg_get_N_insgm;
	part->get_N_nodes = nb_msh3trg_get_N_nodes;
	part->get_N_edges = nb_msh3trg_get_N_edges;
	part->get_N_elems = nb_msh3trg_get_N_elems;
	part->node_get_x = nb_msh3trg_node_get_x;
	part->node_get_y = nb_msh3trg_node_get_y;
	part->edge_get_1n = nb_msh3trg_edge_get_1n;
	part->edge_get_2n = nb_msh3trg_edge_get_2n;
	part->elem_get_x = nb_msh3trg_elem_get_x;
	part->elem_get_y = nb_msh3trg_elem_get_y;
	part->elem_get_area = nb_msh3trg_elem_get_area;
	part->elem_face_get_length = nb_msh3trg_elem_face_get_length;
	part->elem_face_get_normal = nb_msh3trg_elem_face_get_normal;
	part->elem_ngb_get_normal = nb_msh3trg_elem_ngb_get_normal;
	part->elem_get_N_adj = nb_msh3trg_elem_get_N_adj;
	part->elem_get_adj = nb_msh3trg_elem_get_adj;
	part->elem_get_N_ngb = nb_msh3trg_elem_get_N_ngb;
	part->elem_get_ngb = nb_msh3trg_elem_get_ngb;
	part->elem_has_ngb = nb_msh3trg_elem_has_ngb;
	part->get_invtx = nb_msh3trg_get_invtx;
	part->insgm_get_N_nodes = nb_msh3trg_insgm_get_N_nodes;
	part->insgm_get_node = nb_msh3trg_insgm_get_node;
	part->load_elem_graph = nb_msh3trg_load_elem_graph;
	part->load_nodal_graph = nb_msh3trg_load_nodal_graph;
	part->load_interelem_graph = nb_msh3trg_load_interelem_graph;
	part->load_from_mesh = nb_msh3trg_load_from_mesh;
	part->get_enveloping_box = nb_msh3trg_get_enveloping_box;
	part->is_vtx_inside = nb_msh3trg_is_vtx_inside;
	part->build_model = nb_msh3trg_build_model;
	part->build_model_disabled_elems =
		nb_msh3trg_build_model_disabled_elems;
}

static void set_msh3trg_graphics_interface(nb_partition_t *part)
{
	part->graphics.draw_wires = nb_msh3trg_draw_wires;
	part->graphics.draw_boundaries = nb_msh3trg_draw_boundaries;
	part->graphics.fill_elems = nb_msh3trg_fill_elems;
	part->graphics.fill_elems_field_on_nodes =
		nb_msh3trg_fill_elems_field_on_nodes;
	part->graphics.fill_elems_field_on_elems =
		nb_msh3trg_fill_elems_field_on_elems;
	part->graphics.fill_elems_classes = nb_msh3trg_fill_elems_classes;
	part->graphics.fill_nodes = nb_msh3trg_fill_nodes;
	part->graphics.fill_nodes_classes = nb_msh3trg_fill_nodes_classes;
}

static void set_mshquad_interface(nb_partition_t *part)
{
	set_mshquad_main_interface(part);
	set_mshquad_graphics_interface(part);
}

static void set_mshquad_main_interface(nb_partition_t *part)
{
	part->init = nb_mshquad_init;
	part->finish = nb_mshquad_finish;
	part->copy = nb_mshquad_copy;
	part->clear = nb_mshquad_clear;
	part->get_N_invtx = nb_mshquad_get_N_invtx;
	part->get_N_insgm = nb_mshquad_get_N_insgm;
	part->get_N_nodes = nb_mshquad_get_N_nodes;
	part->get_N_edges = nb_mshquad_get_N_edges;
	part->get_N_elems = nb_mshquad_get_N_elems;
	part->node_get_x = nb_mshquad_node_get_x;
	part->node_get_y = nb_mshquad_node_get_y;
	part->edge_get_1n = nb_mshquad_edge_get_1n;
	part->edge_get_2n = nb_mshquad_edge_get_2n;
	part->elem_get_x = nb_mshquad_elem_get_x;
	part->elem_get_y = nb_mshquad_elem_get_y;
	part->elem_get_area = nb_mshquad_elem_get_area;
	part->elem_face_get_length = nb_mshquad_elem_face_get_length;
	part->elem_face_get_normal = nb_mshquad_elem_face_get_normal;
	part->elem_ngb_get_normal = nb_mshquad_elem_ngb_get_normal;
	part->elem_get_N_adj = nb_mshquad_elem_get_N_adj;
	part->elem_get_adj = nb_mshquad_elem_get_adj;
	part->elem_get_N_ngb = nb_mshquad_elem_get_N_ngb;
	part->elem_get_ngb = nb_mshquad_elem_get_ngb;
	part->elem_has_ngb = nb_mshquad_elem_has_ngb;
	part->get_invtx = nb_mshquad_get_invtx;
	part->insgm_get_N_nodes = nb_mshquad_insgm_get_N_nodes;
	part->insgm_get_node = nb_mshquad_insgm_get_node;
	part->load_elem_graph = nb_mshquad_load_elem_graph;
	part->load_nodal_graph = nb_mshquad_load_nodal_graph;
	part->load_interelem_graph = nb_mshquad_load_interelem_graph;
	part->load_from_mesh = nb_mshquad_load_from_mesh;
	part->get_enveloping_box = nb_mshquad_get_enveloping_box;
	part->is_vtx_inside = nb_mshquad_is_vtx_inside;
	part->build_model = nb_mshquad_build_model;
	part->build_model_disabled_elems =
		nb_mshquad_build_model_disabled_elems;
}

static void set_mshquad_graphics_interface(nb_partition_t *part)
{
	part->graphics.draw_wires = nb_mshquad_draw_wires;
	part->graphics.draw_boundaries = nb_mshquad_draw_boundaries;
	part->graphics.fill_elems = nb_mshquad_fill_elems;
	part->graphics.fill_elems_field_on_nodes =
		nb_mshquad_fill_elems_field_on_nodes;
	part->graphics.fill_elems_field_on_elems =
		nb_mshquad_fill_elems_field_on_elems;
	part->graphics.fill_elems_classes = nb_mshquad_fill_elems_classes;
	part->graphics.fill_nodes = nb_mshquad_fill_nodes;
	part->graphics.fill_nodes_classes = nb_mshquad_fill_nodes_classes;
}

static void set_mshpoly_interface(nb_partition_t *part)
{
	set_mshpoly_main_interface(part);
	set_mshpoly_graphics_interface(part);
}

static void set_mshpoly_main_interface(nb_partition_t *part)
{
	part->init = nb_mshpoly_init;
	part->finish = nb_mshpoly_finish;
	part->copy = nb_mshpoly_copy;
	part->clear = nb_mshpoly_clear;
	part->get_N_invtx = nb_mshpoly_get_N_invtx;
	part->get_N_insgm = nb_mshpoly_get_N_insgm;
	part->get_N_nodes = nb_mshpoly_get_N_nodes;
	part->get_N_edges = nb_mshpoly_get_N_edges;
	part->get_N_elems = nb_mshpoly_get_N_elems;
	part->node_get_x = nb_mshpoly_node_get_x;
	part->node_get_y = nb_mshpoly_node_get_y;
	part->edge_get_1n = nb_mshpoly_edge_get_1n;
	part->edge_get_2n = nb_mshpoly_edge_get_2n;
	part->elem_get_x = nb_mshpoly_elem_get_x;
	part->elem_get_y = nb_mshpoly_elem_get_y;
	part->elem_get_area = nb_mshpoly_elem_get_area;
	part->elem_face_get_length = nb_mshpoly_elem_face_get_length;
	part->elem_face_get_normal = nb_mshpoly_elem_face_get_normal;
	part->elem_ngb_get_normal = nb_mshpoly_elem_ngb_get_normal;
	part->elem_get_N_adj = nb_mshpoly_elem_get_N_adj;
	part->elem_get_adj = nb_mshpoly_elem_get_adj;
	part->elem_get_N_ngb = nb_mshpoly_elem_get_N_ngb;
	part->elem_get_ngb = nb_mshpoly_elem_get_ngb;
	part->elem_has_ngb = nb_mshpoly_elem_has_ngb;
	part->get_invtx = nb_mshpoly_get_invtx;
	part->insgm_get_N_nodes = nb_mshpoly_insgm_get_N_nodes;
	part->insgm_get_node = nb_mshpoly_insgm_get_node;
	part->load_elem_graph = nb_mshpoly_load_elem_graph;
	part->load_nodal_graph = nb_mshpoly_load_nodal_graph;
	part->load_interelem_graph = nb_mshpoly_load_interelem_graph;
	part->load_from_mesh = nb_mshpoly_load_from_mesh;
	part->get_enveloping_box = nb_mshpoly_get_enveloping_box;
	part->is_vtx_inside = nb_mshpoly_is_vtx_inside;
	part->build_model = nb_mshpoly_build_model;
	part->build_model_disabled_elems =
		nb_mshpoly_build_model_disabled_elems;
}

static void set_mshpoly_graphics_interface(nb_partition_t *part)
{
	part->graphics.draw_wires = nb_mshpoly_draw_wires;
	part->graphics.draw_boundaries = nb_mshpoly_draw_boundaries;
	part->graphics.fill_elems = nb_mshpoly_fill_elems;
	part->graphics.fill_elems_field_on_nodes =
		nb_mshpoly_fill_elems_field_on_nodes;
	part->graphics.fill_elems_field_on_elems =
		nb_mshpoly_fill_elems_field_on_elems;
	part->graphics.fill_elems_classes = nb_mshpoly_fill_elems_classes;
	part->graphics.fill_nodes = nb_mshpoly_fill_nodes;
	part->graphics.fill_nodes_classes = nb_mshpoly_fill_nodes_classes;
}

static void set_mshpack_interface(nb_partition_t *part)
{
	set_mshpack_main_interface(part);
	set_mshpack_graphics_interface(part);
}

static void set_mshpack_main_interface(nb_partition_t *part)
{
	part->init = nb_mshpack_init;
	part->finish = nb_mshpack_finish;
	part->copy = nb_mshpack_copy;
	part->clear = nb_mshpack_clear;
	part->get_N_invtx = nb_mshpack_get_N_invtx;
	part->get_N_insgm = nb_mshpack_get_N_insgm;
	part->get_N_nodes = nb_mshpack_get_N_nodes;
	part->get_N_edges = nb_mshpack_get_N_edges;
	part->get_N_elems = nb_mshpack_get_N_elems;
	part->node_get_x = nb_mshpack_node_get_x;
	part->node_get_y = nb_mshpack_node_get_y;
	part->edge_get_1n = nb_mshpack_edge_get_1n;
	part->edge_get_2n = nb_mshpack_edge_get_2n;
	part->elem_get_x = nb_mshpack_elem_get_x;
	part->elem_get_y = nb_mshpack_elem_get_y;
	part->elem_get_area = nb_mshpack_elem_get_area;
	part->elem_face_get_length = nb_mshpack_elem_face_get_length;
	part->elem_face_get_normal = nb_mshpack_elem_face_get_normal;
	part->elem_ngb_get_normal = nb_mshpack_elem_ngb_get_normal;
	part->elem_get_N_adj = nb_mshpack_elem_get_N_adj;
	part->elem_get_adj = nb_mshpack_elem_get_adj;
	part->elem_get_N_ngb = nb_mshpack_elem_get_N_ngb;
	part->elem_get_ngb = nb_mshpack_elem_get_ngb;
	part->elem_has_ngb = nb_mshpack_elem_has_ngb;
	part->get_invtx = nb_mshpack_get_invtx;
	part->insgm_get_N_nodes = nb_mshpack_insgm_get_N_nodes;
	part->insgm_get_node = nb_mshpack_insgm_get_node;
	part->load_elem_graph = nb_mshpack_load_elem_graph;
	part->load_nodal_graph = nb_mshpack_load_nodal_graph;
	part->load_interelem_graph = nb_mshpack_load_interelem_graph;
	part->load_from_mesh = nb_mshpack_load_from_mesh;
	part->get_enveloping_box = nb_mshpack_get_enveloping_box;
	part->is_vtx_inside = nb_mshpack_is_vtx_inside;
	part->build_model = nb_mshpack_build_model;
	part->build_model_disabled_elems =
		nb_mshpack_build_model_disabled_elems;
}

static void set_mshpack_graphics_interface(nb_partition_t *part)
{
	part->graphics.draw_wires = nb_mshpack_draw_wires;
	part->graphics.draw_boundaries = nb_mshpack_draw_boundaries;
	part->graphics.fill_elems = nb_mshpack_fill_elems;
	part->graphics.fill_elems_field_on_nodes =
		nb_mshpack_fill_elems_field_on_nodes;
	part->graphics.fill_elems_field_on_elems =
		nb_mshpack_fill_elems_field_on_elems;
	part->graphics.fill_elems_classes = nb_mshpack_fill_elems_classes;
	part->graphics.fill_nodes = nb_mshpack_fill_nodes;
	part->graphics.fill_nodes_classes = nb_mshpack_fill_nodes_classes;
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

double nb_partition_node_get_x(const nb_partition_t *part, uint32_t id)
{
	return part->node_get_x(part->msh, id);
}

double nb_partition_node_get_y(const nb_partition_t *part, uint32_t id)
{
	return part->node_get_y(part->msh, id);
}

double nb_partition_elem_get_area(const nb_partition_t *part, uint32_t id)
{
	return part->elem_get_area(part->msh, id);
}

double nb_partition_elem_face_get_length(const nb_partition_t *part,
					 uint32_t elem_id, uint16_t face_id)
{
	return part->elem_face_get_length(part->msh, elem_id, face_id);
}

void nb_partition_elem_face_get_midpoint(const nb_partition_t *part,
					 uint32_t elem_id, uint16_t face_id,
					 double w, double midpoint[2])
{
	part->elem_face_get_midpoint(part->msh, elem_id, face_id, w, midpoint);
	/* AQUI VOY: Falta agregar en struct y en implementaciones,
	 falta agregar a set_interface. */
}

double nb_partition_elem_face_get_normal(const nb_partition_t *part,
					 uint32_t elem_id, uint16_t face_id,
					 double normal[2])
{
	return part->elem_face_get_normal(part->msh, elem_id,
					  face_id, normal);
}

double nb_partition_elem_ngb_get_normal(const nb_partition_t *part,
					uint32_t elem_id, uint16_t ngb_id,
					double normal[2])
{
	return part->elem_ngb_get_normal(part->msh, elem_id,
					 ngb_id, normal);
}

uint32_t nb_partition_edge_get_1n(const nb_partition_t *part, uint32_t id)
{
	return part->edge_get_1n(part->msh, id);
}

uint32_t nb_partition_edge_get_2n(const nb_partition_t *part, uint32_t id)
{
	return part->edge_get_2n(part->msh, id);
}

double nb_partition_elem_get_x(const nb_partition_t *part, uint32_t id)
{
	return part->elem_get_x(part->msh, id);
}

double nb_partition_elem_get_y(const nb_partition_t *part, uint32_t id)
{
	return part->elem_get_y(part->msh, id);
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

bool nb_partition_elem_has_ngb(const nb_partition_t *part, uint32_t elem_id,
			       uint16_t face_id)
{
	return part->elem_has_ngb(part->msh, elem_id, face_id);
}

uint32_t nb_partition_get_invtx(const nb_partition_t *part, uint32_t id)
{
	return part->get_invtx(part->msh, id);
}

uint32_t nb_partition_insgm_get_N_nodes(const nb_partition_t *part,
					  uint32_t id)
{
	return part->insgm_get_N_nodes(part->msh, id);
}

uint32_t nb_partition_insgm_get_N_subsgm(const nb_partition_t *part,
					 uint32_t id)
{
	uint32_t N = part->insgm_get_N_nodes(part->msh, id);
	return N - 1;
}

uint32_t nb_partition_insgm_get_node(const nb_partition_t *part,
				       uint32_t sgm_id, uint32_t node_id)
{
	return part->insgm_get_node(part->msh, sgm_id, node_id);
}

double nb_partition_insgm_get_length(const nb_partition_t *part,
				     uint32_t sgm_id)
{
	uint32_t last_vtx = part->insgm_get_N_nodes(part->msh, sgm_id) - 1;
	uint32_t n1 = part->insgm_get_node(part->msh, sgm_id, 0);
	uint32_t n2 = part->insgm_get_node(part->msh, sgm_id, last_vtx);
	double x1 = part->node_get_x(part->msh, n1);
	double y1 = part->node_get_y(part->msh, n1);
	double x2 = part->node_get_x(part->msh, n2);
	double y2 = part->node_get_y(part->msh, n2);
	return sqrt(POW2(x1 - x2) + POW2(y1 - y2));
}

double nb_partition_insgm_subsgm_get_length(const nb_partition_t *part,
					    uint32_t sgm_id,
					    uint32_t subsgm_id)
{
	uint32_t n1 = part->insgm_get_node(part->msh, sgm_id, subsgm_id);
	uint32_t n2 = part->insgm_get_node(part->msh, sgm_id, subsgm_id + 1);
	double x1 = part->node_get_x(part->msh, n1);
	double y1 = part->node_get_y(part->msh, n1);
	double x2 = part->node_get_x(part->msh, n2);
	double y2 = part->node_get_y(part->msh, n2);
	return sqrt(POW2(x1 - x2) + POW2(y1 - y2));
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
				 nb_mesh_t *mesh)
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
