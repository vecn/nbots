#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "nb/graph_bot.h"
#include "nb/graphics_bot.h"
#include "nb/geometric_bot/mesh/tessellator2D.h"
#include "nb/geometric_bot/mesh/mesh2D/elements2D/msh3trg.h"
#include "nb/geometric_bot/mesh/mesh2D/elements2D/msh3trg_draw.h"
#include "nb/geometric_bot/mesh/mesh2D/elements2D/mshquad.h"
#include "nb/geometric_bot/mesh/mesh2D/elements2D/mshquad_draw.h"
#include "nb/geometric_bot/mesh/mesh2D/elements2D/mshpoly.h"
#include "nb/geometric_bot/mesh/mesh2D/elements2D/mshpoly_draw.h"
#include "nb/geometric_bot/mesh/mesh2D/elements2D/mshpack.h"
#include "nb/geometric_bot/mesh/mesh2D/elements2D/mshpack_draw.h"

#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/geometric_bot/mesh/mesh2D/info.h"
#include "mesh2D_struct.h"

#define POW2(a) ((a)*(a))

static void set_msh_interface(nb_mesh2D_t *part, nb_mesh2D_type  type);

static void set_msh3trg_interface(nb_mesh2D_t *part);
static void set_msh3trg_main_interface(nb_mesh2D_t *part);
static void set_msh3trg_graphics_interface(nb_mesh2D_t *part);

static void set_mshquad_interface(nb_mesh2D_t *part);
static void set_mshquad_main_interface(nb_mesh2D_t *part);
static void set_mshquad_graphics_interface(nb_mesh2D_t *part);

static void set_mshpoly_interface(nb_mesh2D_t *part);
static void set_mshpoly_main_interface(nb_mesh2D_t *part);
static void set_mshpoly_graphics_interface(nb_mesh2D_t *part);

static void set_mshpack_interface(nb_mesh2D_t *part);
static void set_mshpack_main_interface(nb_mesh2D_t *part);
static void set_mshpack_graphics_interface(nb_mesh2D_t *part);

uint32_t nb_mesh2D_get_memsize(nb_mesh2D_type  type)
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
	return mem + sizeof(nb_mesh2D_t);
}

void nb_mesh2D_init(nb_mesh2D_t *part, nb_mesh2D_type  type)
{
	char *memblock = (void*) part;
	part->msh = (void*) (memblock + sizeof(nb_mesh2D_t));
	part->type = type;
	set_msh_interface(part, type);
	part->init(part->msh);
}

static void set_msh_interface(nb_mesh2D_t *part, nb_mesh2D_type  type)
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

void nb_mesh2D_init_from_msh(nb_mesh2D_t *part, void *msh,
				nb_mesh2D_type  type)
{
	char *memblock = (void*) part;
	part->msh = msh;
	part->type = type;
	set_msh_interface(part, type);
}

static void set_msh3trg_interface(nb_mesh2D_t *part)
{
	set_msh3trg_main_interface(part);
	set_msh3trg_graphics_interface(part);
}

static void set_msh3trg_main_interface(nb_mesh2D_t *part)
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
	part->edge_get_midpoint = nb_msh3trg_edge_get_midpoint;
	part->edge_get_normal = nb_msh3trg_edge_get_normal;
	part->elem_get_x = nb_msh3trg_elem_get_x;
	part->elem_get_y = nb_msh3trg_elem_get_y;
	part->elem_get_area = nb_msh3trg_elem_get_area;
	part->elem_get_radius = nb_msh3trg_elem_get_radius;
	part->elem_get_apotem = nb_msh3trg_elem_get_apotem;
	part->elem_find_edge = nb_msh3trg_elem_find_edge;
	part->elem_face_get_length = nb_msh3trg_elem_face_get_length;
	part->elem_face_get_normal = nb_msh3trg_elem_face_get_normal;
	part->elem_ngb_get_normal = nb_msh3trg_elem_ngb_get_normal;
	part->elem_get_N_adj = nb_msh3trg_elem_get_N_adj;
	part->elem_get_adj = nb_msh3trg_elem_get_adj;
	part->elem_get_ngb = nb_msh3trg_elem_get_ngb;
	part->elem_has_ngb = nb_msh3trg_elem_has_ngb;
	part->get_invtx = nb_msh3trg_get_invtx;
	part->insgm_get_N_nodes = nb_msh3trg_insgm_get_N_nodes;
	part->insgm_get_node = nb_msh3trg_insgm_get_node;
	part->load_from_mesh = nb_msh3trg_load_from_mesh;
	part->set_nodal_permutation = nb_msh3trg_set_nodal_permutation;
	part->get_enveloping_box = nb_msh3trg_get_enveloping_box;
	part->is_vtx_inside = nb_msh3trg_is_vtx_inside;
	part->build_model = nb_msh3trg_build_model;
	part->build_model_disabled_elems =
		nb_msh3trg_build_model_disabled_elems;
}

static void set_msh3trg_graphics_interface(nb_mesh2D_t *part)
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
	part->graphics.draw_level_set = nb_msh3trg_draw_level_set;
}

static void set_mshquad_interface(nb_mesh2D_t *part)
{
	set_mshquad_main_interface(part);
	set_mshquad_graphics_interface(part);
}

static void set_mshquad_main_interface(nb_mesh2D_t *part)
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
	part->edge_get_midpoint = nb_mshquad_edge_get_midpoint;
	part->edge_get_normal = nb_mshquad_edge_get_normal;
	part->elem_get_x = nb_mshquad_elem_get_x;
	part->elem_get_y = nb_mshquad_elem_get_y;
	part->elem_get_area = nb_mshquad_elem_get_area;
	part->elem_get_radius = nb_mshquad_elem_get_radius;
	part->elem_get_apotem = nb_mshquad_elem_get_apotem;
	part->elem_find_edge = nb_mshquad_elem_find_edge;
	part->elem_face_get_length = nb_mshquad_elem_face_get_length;
	part->elem_face_get_normal = nb_mshquad_elem_face_get_normal;
	part->elem_ngb_get_normal = nb_mshquad_elem_ngb_get_normal;
	part->elem_get_N_adj = nb_mshquad_elem_get_N_adj;
	part->elem_get_adj = nb_mshquad_elem_get_adj;
	part->elem_get_ngb = nb_mshquad_elem_get_ngb;
	part->elem_has_ngb = nb_mshquad_elem_has_ngb;
	part->get_invtx = nb_mshquad_get_invtx;
	part->insgm_get_N_nodes = nb_mshquad_insgm_get_N_nodes;
	part->insgm_get_node = nb_mshquad_insgm_get_node;
	part->load_from_mesh = nb_mshquad_load_from_mesh;
	part->set_nodal_permutation = nb_mshquad_set_nodal_permutation;
	part->get_enveloping_box = nb_mshquad_get_enveloping_box;
	part->is_vtx_inside = nb_mshquad_is_vtx_inside;
	part->build_model = nb_mshquad_build_model;
	part->build_model_disabled_elems =
		nb_mshquad_build_model_disabled_elems;
}

static void set_mshquad_graphics_interface(nb_mesh2D_t *part)
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
	part->graphics.draw_level_set = nb_mshquad_draw_level_set;
}

static void set_mshpoly_interface(nb_mesh2D_t *part)
{
	set_mshpoly_main_interface(part);
	set_mshpoly_graphics_interface(part);
}

static void set_mshpoly_main_interface(nb_mesh2D_t *part)
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
	part->edge_get_midpoint = nb_mshpoly_edge_get_midpoint;
	part->edge_get_normal = nb_mshpoly_edge_get_normal;
	part->elem_get_x = nb_mshpoly_elem_get_x;
	part->elem_get_y = nb_mshpoly_elem_get_y;
	part->elem_get_area = nb_mshpoly_elem_get_area;
	part->elem_get_radius = nb_mshpoly_elem_get_radius;
	part->elem_get_apotem = nb_mshpoly_elem_get_apotem;
	part->elem_find_edge = nb_mshpoly_elem_find_edge;
	part->elem_face_get_length = nb_mshpoly_elem_face_get_length;
	part->elem_face_get_normal = nb_mshpoly_elem_face_get_normal;
	part->elem_ngb_get_normal = nb_mshpoly_elem_ngb_get_normal;
	part->elem_get_N_adj = nb_mshpoly_elem_get_N_adj;
	part->elem_get_adj = nb_mshpoly_elem_get_adj;
	part->elem_get_ngb = nb_mshpoly_elem_get_ngb;
	part->elem_has_ngb = nb_mshpoly_elem_has_ngb;
	part->get_invtx = nb_mshpoly_get_invtx;
	part->insgm_get_N_nodes = nb_mshpoly_insgm_get_N_nodes;
	part->insgm_get_node = nb_mshpoly_insgm_get_node;
	part->load_from_mesh = nb_mshpoly_load_from_mesh;
	part->set_nodal_permutation = nb_mshpoly_set_nodal_permutation;
	part->get_enveloping_box = nb_mshpoly_get_enveloping_box;
	part->is_vtx_inside = nb_mshpoly_is_vtx_inside;
	part->build_model = nb_mshpoly_build_model;
	part->build_model_disabled_elems =
		nb_mshpoly_build_model_disabled_elems;
}

static void set_mshpoly_graphics_interface(nb_mesh2D_t *part)
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
	part->graphics.draw_level_set = nb_mshpoly_draw_level_set;
}

static void set_mshpack_interface(nb_mesh2D_t *part)
{
	set_mshpack_main_interface(part);
	set_mshpack_graphics_interface(part);
}

static void set_mshpack_main_interface(nb_mesh2D_t *part)
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
	part->edge_get_midpoint = nb_mshpack_edge_get_midpoint;
	part->edge_get_normal = nb_mshpack_edge_get_normal;
	part->elem_get_x = nb_mshpack_elem_get_x;
	part->elem_get_y = nb_mshpack_elem_get_y;
	part->elem_get_area = nb_mshpack_elem_get_area;
	part->elem_get_radius = nb_mshpack_elem_get_radius;
	part->elem_get_apotem = nb_mshpack_elem_get_apotem;
	part->elem_find_edge = nb_mshpack_elem_find_edge;
	part->elem_face_get_length = nb_mshpack_elem_face_get_length;
	part->elem_face_get_normal = nb_mshpack_elem_face_get_normal;
	part->elem_ngb_get_normal = nb_mshpack_elem_ngb_get_normal;
	part->elem_get_N_adj = nb_mshpack_elem_get_N_adj;
	part->elem_get_adj = nb_mshpack_elem_get_adj;
	part->elem_get_ngb = nb_mshpack_elem_get_ngb;
	part->elem_has_ngb = nb_mshpack_elem_has_ngb;
	part->get_invtx = nb_mshpack_get_invtx;
	part->insgm_get_N_nodes = nb_mshpack_insgm_get_N_nodes;
	part->insgm_get_node = nb_mshpack_insgm_get_node;
	part->load_from_mesh = nb_mshpack_load_from_mesh;
	part->set_nodal_permutation = nb_mshpack_set_nodal_permutation;
	part->get_enveloping_box = nb_mshpack_get_enveloping_box;
	part->is_vtx_inside = nb_mshpack_is_vtx_inside;
	part->build_model = nb_mshpack_build_model;
	part->build_model_disabled_elems =
		nb_mshpack_build_model_disabled_elems;
}

static void set_mshpack_graphics_interface(nb_mesh2D_t *part)
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
	part->graphics.draw_level_set = nb_mshpack_draw_level_set;
}

void nb_mesh2D_copy(nb_mesh2D_t *part, const nb_mesh2D_t* srcpart)
{
	nb_mesh2D_init(part, srcpart->type);
	part->copy(part->msh, srcpart->msh);
}

void nb_mesh2D_finish(nb_mesh2D_t *part)
{
	part->finish(part->msh);
}

nb_mesh2D_t* nb_mesh2D_create(nb_mesh2D_type type)
{
	uint32_t memsize = nb_mesh2D_get_memsize(type);
	nb_mesh2D_t *part = nb_allocate_mem(memsize);
	nb_mesh2D_init(part, type);
	return part;
}

nb_mesh2D_t* nb_mesh2D_clone(nb_mesh2D_t* part)
{
	uint32_t memsize = nb_mesh2D_get_memsize(part->type);
	nb_mesh2D_t *clone = nb_allocate_mem(memsize);
	nb_mesh2D_copy(clone, part);
	return clone;
}

void nb_mesh2D_clear(nb_mesh2D_t* part)
{
	part->clear(part->msh);
}

void nb_mesh2D_destroy(nb_mesh2D_t* part)
{
	nb_mesh2D_finish(part);
	nb_free_mem(part);
}

nb_mesh2D_type nb_mesh2D_get_type(const nb_mesh2D_t *part)
{
	return part->type;
}

uint32_t nb_mesh2D_get_N_invtx(const nb_mesh2D_t *part)
{
	return part->get_N_invtx(part->msh);
}

uint32_t nb_mesh2D_get_N_insgm(const nb_mesh2D_t *part)
{
	return part->get_N_insgm(part->msh);
}

uint32_t nb_mesh2D_get_N_nodes(const nb_mesh2D_t *part)
{
	return part->get_N_nodes(part->msh);
}

uint32_t nb_mesh2D_get_N_edges(const nb_mesh2D_t *part)
{
	return part->get_N_edges(part->msh);
}

uint32_t nb_mesh2D_get_N_elems(const nb_mesh2D_t *part)
{
	return part->get_N_elems(part->msh);
}

double nb_mesh2D_node_get_x(const nb_mesh2D_t *part, uint32_t id)
{
	return part->node_get_x(part->msh, id);
}

double nb_mesh2D_node_get_y(const nb_mesh2D_t *part, uint32_t id)
{
	return part->node_get_y(part->msh, id);
}

uint32_t nb_mesh2D_edge_get_1n(const nb_mesh2D_t *part, uint32_t id)
{
	return part->edge_get_1n(part->msh, id);
}

uint32_t nb_mesh2D_edge_get_2n(const nb_mesh2D_t *part, uint32_t id)
{
	return part->edge_get_2n(part->msh, id);
}

void nb_mesh2D_edge_get_midpoint(const nb_mesh2D_t *part,
				    uint32_t face_id, double w,
				    double midpoint[2])
{
	part->edge_get_midpoint(part->msh, face_id, w, midpoint);
}

double nb_mesh2D_edge_get_normal(const nb_mesh2D_t *part,
				    uint32_t face_id, double normal[2])
{
	return part->edge_get_normal(part->msh, face_id, normal);
}

double nb_mesh2D_edge_get_length(const nb_mesh2D_t *part,
				    uint32_t face_id)
{
	double nf[2];
	double length = part->edge_get_normal(part->msh, face_id, nf);
	return length;
}

double nb_mesh2D_elem_get_area(const nb_mesh2D_t *part, uint32_t id)
{
	return part->elem_get_area(part->msh, id);
}

double nb_mesh2D_elem_get_radius(const nb_mesh2D_t *part, uint32_t id)
{
	return part->elem_get_radius(part->msh, id);
}

double nb_mesh2D_elem_get_apotem(const nb_mesh2D_t *part, uint32_t id)
{
	return part->elem_get_apotem(part->msh, id);
}

uint32_t nb_mesh2D_elem_find_edge(const nb_mesh2D_t *part, uint32_t id,
				     uint16_t local_face_id)
{
	return part->elem_find_edge(part->msh, id, local_face_id);
}

double nb_mesh2D_elem_face_get_length(const nb_mesh2D_t *part,
					 uint32_t elem_id, uint16_t face_id)
{
	return part->elem_face_get_length(part->msh, elem_id, face_id);
}

double nb_mesh2D_elem_face_get_normal(const nb_mesh2D_t *part,
					 uint32_t elem_id, uint16_t face_id,
					 double normal[2])
{
	return part->elem_face_get_normal(part->msh, elem_id,
					  face_id, normal);
}

uint32_t nb_mesh2D_elem_face_get_left_ngb(const nb_mesh2D_t *part,
					     uint32_t elem_id,
					     uint16_t face_id)
{
	uint16_t N_adj = nb_mesh2D_elem_get_N_adj(part, elem_id);

	uint32_t left_elem;
	if (face_id < N_adj) {
		face_id = (face_id + 1) % N_adj;
		left_elem = nb_mesh2D_elem_get_ngb(part, elem_id,
						      face_id);
	} else {
		left_elem = nb_mesh2D_get_N_elems(part);
	}
	return left_elem;	
}

uint32_t nb_mesh2D_elem_face_get_right_ngb(const nb_mesh2D_t *part,
					      uint32_t elem_id,
					      uint16_t face_id)
{
	uint16_t N_adj = nb_mesh2D_elem_get_N_adj(part, elem_id);

	uint32_t right_elem;
	if (face_id < N_adj) {
		face_id = (face_id == 0)?(N_adj-1):(face_id-1);
		right_elem = nb_mesh2D_elem_get_ngb(part, elem_id,
						       face_id);
	} else {
		right_elem = nb_mesh2D_get_N_elems(part);
	}
	return right_elem;	
}

double nb_mesh2D_elem_ngb_get_normal(const nb_mesh2D_t *part,
					uint32_t elem_id, uint16_t ngb_id,
					double normal[2])
{
	return part->elem_ngb_get_normal(part->msh, elem_id,
					 ngb_id, normal);
}

uint16_t nb_mesh2D_elem_ngb_get_face(const nb_mesh2D_t *part,
					uint32_t elem_id, uint32_t ngb_id)
{
	uint16_t N_adj = nb_mesh2D_elem_get_N_adj(part, elem_id);
	uint16_t face_id = N_adj;
	for (uint16_t i = 0; i < N_adj; i++) {
		uint32_t ingb = nb_mesh2D_elem_get_ngb(part, elem_id, i);
		if (ingb == ngb_id) {
			face_id = i;
			break;
		}
	}
	return face_id;
}

double nb_mesh2D_elem_get_x(const nb_mesh2D_t *part, uint32_t id)
{
	return part->elem_get_x(part->msh, id);
}

double nb_mesh2D_elem_get_y(const nb_mesh2D_t *part, uint32_t id)
{
	return part->elem_get_y(part->msh, id);
}

uint32_t nb_mesh2D_elem_get_N_adj(const nb_mesh2D_t *part, uint32_t id)
{
	return part->elem_get_N_adj(part->msh, id);
}

uint32_t nb_mesh2D_elem_get_adj(const nb_mesh2D_t *part,
				   uint32_t elem_id, uint8_t adj_id)
{
	return part->elem_get_adj(part->msh, elem_id, adj_id);
}

uint32_t nb_mesh2D_elem_get_ngb(const nb_mesh2D_t *part,
				   uint32_t elem_id, uint8_t ngb_id)
{
	return part->elem_get_ngb(part->msh, elem_id, ngb_id);
}

bool nb_mesh2D_elem_has_ngb(const nb_mesh2D_t *part, uint32_t elem_id,
			       uint16_t face_id)
{
	return part->elem_has_ngb(part->msh, elem_id, face_id);
}

uint32_t nb_mesh2D_get_invtx(const nb_mesh2D_t *part, uint32_t id)
{
	return part->get_invtx(part->msh, id);
}

uint32_t nb_mesh2D_insgm_get_N_nodes(const nb_mesh2D_t *part,
					  uint32_t id)
{
	return part->insgm_get_N_nodes(part->msh, id);
}

uint32_t nb_mesh2D_insgm_get_N_subsgm(const nb_mesh2D_t *part,
					 uint32_t id)
{
	uint32_t N = part->insgm_get_N_nodes(part->msh, id);
	return N - 1;
}

uint32_t nb_mesh2D_insgm_get_node(const nb_mesh2D_t *part,
				       uint32_t sgm_id, uint32_t node_id)
{
	return part->insgm_get_node(part->msh, sgm_id, node_id);
}

double nb_mesh2D_insgm_get_length(const nb_mesh2D_t *part,
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

double nb_mesh2D_insgm_subsgm_get_length(const nb_mesh2D_t *part,
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

void nb_mesh2D_load_from_mesh(nb_mesh2D_t *part,
				 nb_tessellator2D_t *mesh)
{
	part->load_from_mesh(part->msh, mesh);
}

void nb_mesh2D_set_nodal_permutation(nb_mesh2D_t *part,
					const uint32_t *perm)
{
	part->set_nodal_permutation(part->msh, perm);
}

void nb_mesh2D_get_enveloping_box(const nb_mesh2D_t *part,
				     double box[4])
{
	part->get_enveloping_box(part->msh, box);
}

bool nb_mesh2D_is_vtx_inside(const nb_mesh2D_t *part,
				double x, double y)
{
	return part->is_vtx_inside(part->msh, x, y);
}

void nb_mesh2D_build_model(const nb_mesh2D_t *part, nb_model_t *model)
{
	part->build_model(part->msh, model);
}

void nb_mesh2D_build_model_disabled_elems
			(const nb_mesh2D_t *part,
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
