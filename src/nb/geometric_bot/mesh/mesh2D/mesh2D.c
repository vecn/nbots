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

static void set_msh_interface(nb_mesh2D_t *mesh, nb_mesh2D_type  type);

static void set_msh3trg_interface(nb_mesh2D_t *mesh);
static void set_msh3trg_main_interface(nb_mesh2D_t *mesh);
static void set_msh3trg_graphics_interface(nb_mesh2D_t *mesh);

static void set_mshquad_interface(nb_mesh2D_t *mesh);
static void set_mshquad_main_interface(nb_mesh2D_t *mesh);
static void set_mshquad_graphics_interface(nb_mesh2D_t *mesh);

static void set_mshpoly_interface(nb_mesh2D_t *mesh);
static void set_mshpoly_main_interface(nb_mesh2D_t *mesh);
static void set_mshpoly_graphics_interface(nb_mesh2D_t *mesh);

static void set_mshpack_interface(nb_mesh2D_t *mesh);
static void set_mshpack_main_interface(nb_mesh2D_t *mesh);
static void set_mshpack_graphics_interface(nb_mesh2D_t *mesh);

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

void nb_mesh2D_init(nb_mesh2D_t *mesh, nb_mesh2D_type  type)
{
	char *memblock = (void*) mesh;
	mesh->msh = (void*) (memblock + sizeof(nb_mesh2D_t));
	mesh->type = type;
	set_msh_interface(mesh, type);
	mesh->init(mesh->msh);
}

static void set_msh_interface(nb_mesh2D_t *mesh, nb_mesh2D_type  type)
{
	switch (type) {
	case NB_TRIAN:
		set_msh3trg_interface(mesh);
		break;
	case NB_QUAD:
		set_mshquad_interface(mesh);
		break;
	case NB_POLY:
		set_mshpoly_interface(mesh);
		break;
	case NB_DISK:
		set_mshpack_interface(mesh);
		break;
	default:
		set_msh3trg_interface(mesh);
		break;
	}
}

void nb_mesh2D_init_from_msh(nb_mesh2D_t *mesh, void *msh,
				nb_mesh2D_type  type)
{
	char *memblock = (void*) mesh;
	mesh->msh = msh;
	mesh->type = type;
	set_msh_interface(mesh, type);
}

static void set_msh3trg_interface(nb_mesh2D_t *mesh)
{
	set_msh3trg_main_interface(mesh);
	set_msh3trg_graphics_interface(mesh);
}

static void set_msh3trg_main_interface(nb_mesh2D_t *mesh)
{
	mesh->init = nb_msh3trg_init;
	mesh->finish = nb_msh3trg_finish;
	mesh->copy = nb_msh3trg_copy;
	mesh->clear = nb_msh3trg_clear;
	mesh->get_N_invtx = nb_msh3trg_get_N_invtx;
	mesh->get_N_insgm = nb_msh3trg_get_N_insgm;
	mesh->get_N_nodes = nb_msh3trg_get_N_nodes;
	mesh->get_N_edges = nb_msh3trg_get_N_edges;
	mesh->get_N_elems = nb_msh3trg_get_N_elems;
	mesh->node_get_x = nb_msh3trg_node_get_x;
	mesh->node_get_y = nb_msh3trg_node_get_y;
	mesh->edge_get_1n = nb_msh3trg_edge_get_1n;
	mesh->edge_get_2n = nb_msh3trg_edge_get_2n;
	mesh->edge_get_midpoint = nb_msh3trg_edge_get_midpoint;
	mesh->edge_get_normal = nb_msh3trg_edge_get_normal;
	mesh->elem_get_x = nb_msh3trg_elem_get_x;
	mesh->elem_get_y = nb_msh3trg_elem_get_y;
	mesh->elem_get_area = nb_msh3trg_elem_get_area;
	mesh->elem_get_radius = nb_msh3trg_elem_get_radius;
	mesh->elem_get_apotem = nb_msh3trg_elem_get_apotem;
	mesh->elem_find_edge = nb_msh3trg_elem_find_edge;
	mesh->elem_face_get_length = nb_msh3trg_elem_face_get_length;
	mesh->elem_face_get_normal = nb_msh3trg_elem_face_get_normal;
	mesh->elem_ngb_get_normal = nb_msh3trg_elem_ngb_get_normal;
	mesh->elem_get_N_adj = nb_msh3trg_elem_get_N_adj;
	mesh->elem_get_adj = nb_msh3trg_elem_get_adj;
	mesh->elem_get_ngb = nb_msh3trg_elem_get_ngb;
	mesh->elem_has_ngb = nb_msh3trg_elem_has_ngb;
	mesh->elem_is_boundary = nb_msh3trg_elem_is_boundary;
	mesh->get_invtx = nb_msh3trg_get_invtx;
	mesh->insgm_get_N_nodes = nb_msh3trg_insgm_get_N_nodes;
	mesh->insgm_get_node = nb_msh3trg_insgm_get_node;
	mesh->load_from_mesh = nb_msh3trg_load_from_tessellator2D;
	mesh->set_nodal_permutation = nb_msh3trg_set_nodal_permutation;
	mesh->get_enveloping_box = nb_msh3trg_get_enveloping_box;
	mesh->is_vtx_inside = nb_msh3trg_is_vtx_inside;
	mesh->build_model = nb_msh3trg_build_model;
	mesh->build_model_disabled_elems =
		nb_msh3trg_build_model_disabled_elems;
}

static void set_msh3trg_graphics_interface(nb_mesh2D_t *mesh)
{
	mesh->graphics.draw_wires = nb_msh3trg_draw_wires;
	mesh->graphics.draw_boundaries = nb_msh3trg_draw_boundaries;
	mesh->graphics.fill_elems = nb_msh3trg_fill_elems;
	mesh->graphics.fill_elems_field_on_nodes =
		nb_msh3trg_fill_elems_field_on_nodes;
	mesh->graphics.fill_elems_field_on_elems =
		nb_msh3trg_fill_elems_field_on_elems;
	mesh->graphics.fill_elems_classes = nb_msh3trg_fill_elems_classes;
	mesh->graphics.fill_nodes = nb_msh3trg_fill_nodes;
	mesh->graphics.fill_nodes_classes = nb_msh3trg_fill_nodes_classes;
	mesh->graphics.draw_level_set = nb_msh3trg_draw_level_set;
}

static void set_mshquad_interface(nb_mesh2D_t *mesh)
{
	set_mshquad_main_interface(mesh);
	set_mshquad_graphics_interface(mesh);
}

static void set_mshquad_main_interface(nb_mesh2D_t *mesh)
{
	mesh->init = nb_mshquad_init;
	mesh->finish = nb_mshquad_finish;
	mesh->copy = nb_mshquad_copy;
	mesh->clear = nb_mshquad_clear;
	mesh->get_N_invtx = nb_mshquad_get_N_invtx;
	mesh->get_N_insgm = nb_mshquad_get_N_insgm;
	mesh->get_N_nodes = nb_mshquad_get_N_nodes;
	mesh->get_N_edges = nb_mshquad_get_N_edges;
	mesh->get_N_elems = nb_mshquad_get_N_elems;
	mesh->node_get_x = nb_mshquad_node_get_x;
	mesh->node_get_y = nb_mshquad_node_get_y;
	mesh->edge_get_1n = nb_mshquad_edge_get_1n;
	mesh->edge_get_2n = nb_mshquad_edge_get_2n;
	mesh->edge_get_midpoint = nb_mshquad_edge_get_midpoint;
	mesh->edge_get_normal = nb_mshquad_edge_get_normal;
	mesh->elem_get_x = nb_mshquad_elem_get_x;
	mesh->elem_get_y = nb_mshquad_elem_get_y;
	mesh->elem_get_area = nb_mshquad_elem_get_area;
	mesh->elem_get_radius = nb_mshquad_elem_get_radius;
	mesh->elem_get_apotem = nb_mshquad_elem_get_apotem;
	mesh->elem_find_edge = nb_mshquad_elem_find_edge;
	mesh->elem_face_get_length = nb_mshquad_elem_face_get_length;
	mesh->elem_face_get_normal = nb_mshquad_elem_face_get_normal;
	mesh->elem_ngb_get_normal = nb_mshquad_elem_ngb_get_normal;
	mesh->elem_get_N_adj = nb_mshquad_elem_get_N_adj;
	mesh->elem_get_adj = nb_mshquad_elem_get_adj;
	mesh->elem_get_ngb = nb_mshquad_elem_get_ngb;
	mesh->elem_has_ngb = nb_mshquad_elem_has_ngb;
	mesh->elem_is_boundary = nb_mshquad_elem_is_boundary;
	mesh->get_invtx = nb_mshquad_get_invtx;
	mesh->insgm_get_N_nodes = nb_mshquad_insgm_get_N_nodes;
	mesh->insgm_get_node = nb_mshquad_insgm_get_node;
	mesh->load_from_mesh = nb_mshquad_load_from_tessellator2D;
	mesh->set_nodal_permutation = nb_mshquad_set_nodal_permutation;
	mesh->get_enveloping_box = nb_mshquad_get_enveloping_box;
	mesh->is_vtx_inside = nb_mshquad_is_vtx_inside;
	mesh->build_model = nb_mshquad_build_model;
	mesh->build_model_disabled_elems =
		nb_mshquad_build_model_disabled_elems;
}

static void set_mshquad_graphics_interface(nb_mesh2D_t *mesh)
{
	mesh->graphics.draw_wires = nb_mshquad_draw_wires;
	mesh->graphics.draw_boundaries = nb_mshquad_draw_boundaries;
	mesh->graphics.fill_elems = nb_mshquad_fill_elems;
	mesh->graphics.fill_elems_field_on_nodes =
		nb_mshquad_fill_elems_field_on_nodes;
	mesh->graphics.fill_elems_field_on_elems =
		nb_mshquad_fill_elems_field_on_elems;
	mesh->graphics.fill_elems_classes = nb_mshquad_fill_elems_classes;
	mesh->graphics.fill_nodes = nb_mshquad_fill_nodes;
	mesh->graphics.fill_nodes_classes = nb_mshquad_fill_nodes_classes;
	mesh->graphics.draw_level_set = nb_mshquad_draw_level_set;
}

static void set_mshpoly_interface(nb_mesh2D_t *mesh)
{
	set_mshpoly_main_interface(mesh);
	set_mshpoly_graphics_interface(mesh);
}

static void set_mshpoly_main_interface(nb_mesh2D_t *mesh)
{
	mesh->init = nb_mshpoly_init;
	mesh->finish = nb_mshpoly_finish;
	mesh->copy = nb_mshpoly_copy;
	mesh->clear = nb_mshpoly_clear;
	mesh->get_N_invtx = nb_mshpoly_get_N_invtx;
	mesh->get_N_insgm = nb_mshpoly_get_N_insgm;
	mesh->get_N_nodes = nb_mshpoly_get_N_nodes;
	mesh->get_N_edges = nb_mshpoly_get_N_edges;
	mesh->get_N_elems = nb_mshpoly_get_N_elems;
	mesh->node_get_x = nb_mshpoly_node_get_x;
	mesh->node_get_y = nb_mshpoly_node_get_y;
	mesh->edge_get_1n = nb_mshpoly_edge_get_1n;
	mesh->edge_get_2n = nb_mshpoly_edge_get_2n;
	mesh->edge_get_midpoint = nb_mshpoly_edge_get_midpoint;
	mesh->edge_get_normal = nb_mshpoly_edge_get_normal;
	mesh->elem_get_x = nb_mshpoly_elem_get_x;
	mesh->elem_get_y = nb_mshpoly_elem_get_y;
	mesh->elem_get_area = nb_mshpoly_elem_get_area;
	mesh->elem_get_radius = nb_mshpoly_elem_get_radius;
	mesh->elem_get_apotem = nb_mshpoly_elem_get_apotem;
	mesh->elem_find_edge = nb_mshpoly_elem_find_edge;
	mesh->elem_face_get_length = nb_mshpoly_elem_face_get_length;
	mesh->elem_face_get_normal = nb_mshpoly_elem_face_get_normal;
	mesh->elem_ngb_get_normal = nb_mshpoly_elem_ngb_get_normal;
	mesh->elem_get_N_adj = nb_mshpoly_elem_get_N_adj;
	mesh->elem_get_adj = nb_mshpoly_elem_get_adj;
	mesh->elem_get_ngb = nb_mshpoly_elem_get_ngb;
	mesh->elem_has_ngb = nb_mshpoly_elem_has_ngb;
	mesh->elem_is_boundary = nb_mshpoly_elem_is_boundary;
	mesh->get_invtx = nb_mshpoly_get_invtx;
	mesh->insgm_get_N_nodes = nb_mshpoly_insgm_get_N_nodes;
	mesh->insgm_get_node = nb_mshpoly_insgm_get_node;
	mesh->load_from_mesh = nb_mshpoly_load_from_tessellator2D;
	mesh->set_nodal_permutation = nb_mshpoly_set_nodal_permutation;
	mesh->get_enveloping_box = nb_mshpoly_get_enveloping_box;
	mesh->is_vtx_inside = nb_mshpoly_is_vtx_inside;
	mesh->build_model = nb_mshpoly_build_model;
	mesh->build_model_disabled_elems =
		nb_mshpoly_build_model_disabled_elems;
}

static void set_mshpoly_graphics_interface(nb_mesh2D_t *mesh)
{
	mesh->graphics.draw_wires = nb_mshpoly_draw_wires;
	mesh->graphics.draw_boundaries = nb_mshpoly_draw_boundaries;
	mesh->graphics.fill_elems = nb_mshpoly_fill_elems;
	mesh->graphics.fill_elems_field_on_nodes =
		nb_mshpoly_fill_elems_field_on_nodes;
	mesh->graphics.fill_elems_field_on_elems =
		nb_mshpoly_fill_elems_field_on_elems;
	mesh->graphics.fill_elems_classes = nb_mshpoly_fill_elems_classes;
	mesh->graphics.fill_nodes = nb_mshpoly_fill_nodes;
	mesh->graphics.fill_nodes_classes = nb_mshpoly_fill_nodes_classes;
	mesh->graphics.draw_level_set = nb_mshpoly_draw_level_set;
}

static void set_mshpack_interface(nb_mesh2D_t *mesh)
{
	set_mshpack_main_interface(mesh);
	set_mshpack_graphics_interface(mesh);
}

static void set_mshpack_main_interface(nb_mesh2D_t *mesh)
{
	mesh->init = nb_mshpack_init;
	mesh->finish = nb_mshpack_finish;
	mesh->copy = nb_mshpack_copy;
	mesh->clear = nb_mshpack_clear;
	mesh->get_N_invtx = nb_mshpack_get_N_invtx;
	mesh->get_N_insgm = nb_mshpack_get_N_insgm;
	mesh->get_N_nodes = nb_mshpack_get_N_nodes;
	mesh->get_N_edges = nb_mshpack_get_N_edges;
	mesh->get_N_elems = nb_mshpack_get_N_elems;
	mesh->node_get_x = nb_mshpack_node_get_x;
	mesh->node_get_y = nb_mshpack_node_get_y;
	mesh->edge_get_1n = nb_mshpack_edge_get_1n;
	mesh->edge_get_2n = nb_mshpack_edge_get_2n;
	mesh->edge_get_midpoint = nb_mshpack_edge_get_midpoint;
	mesh->edge_get_normal = nb_mshpack_edge_get_normal;
	mesh->elem_get_x = nb_mshpack_elem_get_x;
	mesh->elem_get_y = nb_mshpack_elem_get_y;
	mesh->elem_get_area = nb_mshpack_elem_get_area;
	mesh->elem_get_radius = nb_mshpack_elem_get_radius;
	mesh->elem_get_apotem = nb_mshpack_elem_get_apotem;
	mesh->elem_find_edge = nb_mshpack_elem_find_edge;
	mesh->elem_face_get_length = nb_mshpack_elem_face_get_length;
	mesh->elem_face_get_normal = nb_mshpack_elem_face_get_normal;
	mesh->elem_ngb_get_normal = nb_mshpack_elem_ngb_get_normal;
	mesh->elem_get_N_adj = nb_mshpack_elem_get_N_adj;
	mesh->elem_get_adj = nb_mshpack_elem_get_adj;
	mesh->elem_get_ngb = nb_mshpack_elem_get_ngb;
	mesh->elem_has_ngb = nb_mshpack_elem_has_ngb;
	mesh->elem_is_boundary = nb_mshpack_elem_is_boundary;
	mesh->get_invtx = nb_mshpack_get_invtx;
	mesh->insgm_get_N_nodes = nb_mshpack_insgm_get_N_nodes;
	mesh->insgm_get_node = nb_mshpack_insgm_get_node;
	mesh->load_from_mesh = nb_mshpack_load_from_tessellator2D;
	mesh->set_nodal_permutation = nb_mshpack_set_nodal_permutation;
	mesh->get_enveloping_box = nb_mshpack_get_enveloping_box;
	mesh->is_vtx_inside = nb_mshpack_is_vtx_inside;
	mesh->build_model = nb_mshpack_build_model;
	mesh->build_model_disabled_elems =
		nb_mshpack_build_model_disabled_elems;
}

static void set_mshpack_graphics_interface(nb_mesh2D_t *mesh)
{
	mesh->graphics.draw_wires = nb_mshpack_draw_wires;
	mesh->graphics.draw_boundaries = nb_mshpack_draw_boundaries;
	mesh->graphics.fill_elems = nb_mshpack_fill_elems;
	mesh->graphics.fill_elems_field_on_nodes =
		nb_mshpack_fill_elems_field_on_nodes;
	mesh->graphics.fill_elems_field_on_elems =
		nb_mshpack_fill_elems_field_on_elems;
	mesh->graphics.fill_elems_classes = nb_mshpack_fill_elems_classes;
	mesh->graphics.fill_nodes = nb_mshpack_fill_nodes;
	mesh->graphics.fill_nodes_classes = nb_mshpack_fill_nodes_classes;
	mesh->graphics.draw_level_set = nb_mshpack_draw_level_set;
}

void nb_mesh2D_copy(nb_mesh2D_t *mesh, const nb_mesh2D_t* srcmesh)
{
	nb_mesh2D_init(mesh, srcmesh->type);
	mesh->copy(mesh->msh, srcmesh->msh);
}

void nb_mesh2D_finish(nb_mesh2D_t *mesh)
{
	mesh->finish(mesh->msh);
}

nb_mesh2D_t* nb_mesh2D_create(nb_mesh2D_type type)
{
	uint32_t memsize = nb_mesh2D_get_memsize(type);
	nb_mesh2D_t *mesh = nb_allocate_mem(memsize);
	nb_mesh2D_init(mesh, type);
	return mesh;
}

nb_mesh2D_t* nb_mesh2D_clone(nb_mesh2D_t* mesh)
{
	uint32_t memsize = nb_mesh2D_get_memsize(mesh->type);
	nb_mesh2D_t *clone = nb_allocate_mem(memsize);
	nb_mesh2D_copy(clone, mesh);
	return clone;
}

void nb_mesh2D_clear(nb_mesh2D_t* mesh)
{
	mesh->clear(mesh->msh);
}

void nb_mesh2D_destroy(nb_mesh2D_t* mesh)
{
	nb_mesh2D_finish(mesh);
	nb_free_mem(mesh);
}

nb_mesh2D_type nb_mesh2D_get_type(const nb_mesh2D_t *mesh)
{
	return mesh->type;
}

uint32_t nb_mesh2D_get_N_invtx(const nb_mesh2D_t *mesh)
{
	return mesh->get_N_invtx(mesh->msh);
}

uint32_t nb_mesh2D_get_N_insgm(const nb_mesh2D_t *mesh)
{
	return mesh->get_N_insgm(mesh->msh);
}

uint32_t nb_mesh2D_get_N_nodes(const nb_mesh2D_t *mesh)
{
	return mesh->get_N_nodes(mesh->msh);
}

uint32_t nb_mesh2D_get_N_edges(const nb_mesh2D_t *mesh)
{
	return mesh->get_N_edges(mesh->msh);
}

uint32_t nb_mesh2D_get_N_elems(const nb_mesh2D_t *mesh)
{
	return mesh->get_N_elems(mesh->msh);
}

double nb_mesh2D_node_get_x(const nb_mesh2D_t *mesh, uint32_t id)
{
	return mesh->node_get_x(mesh->msh, id);
}

double nb_mesh2D_node_get_y(const nb_mesh2D_t *mesh, uint32_t id)
{
	return mesh->node_get_y(mesh->msh, id);
}

uint32_t nb_mesh2D_edge_get_1n(const nb_mesh2D_t *mesh, uint32_t id)
{
	return mesh->edge_get_1n(mesh->msh, id);
}

uint32_t nb_mesh2D_edge_get_2n(const nb_mesh2D_t *mesh, uint32_t id)
{
	return mesh->edge_get_2n(mesh->msh, id);
}

void nb_mesh2D_edge_get_midpoint(const nb_mesh2D_t *mesh,
				    uint32_t face_id, double w,
				    double midpoint[2])
{
	mesh->edge_get_midpoint(mesh->msh, face_id, w, midpoint);
}

double nb_mesh2D_edge_get_normal(const nb_mesh2D_t *mesh,
				    uint32_t face_id, double normal[2])
{
	return mesh->edge_get_normal(mesh->msh, face_id, normal);
}

double nb_mesh2D_edge_get_length(const nb_mesh2D_t *mesh,
				    uint32_t face_id)
{
	double nf[2];
	double length = mesh->edge_get_normal(mesh->msh, face_id, nf);
	return length;
}

double nb_mesh2D_elem_get_area(const nb_mesh2D_t *mesh, uint32_t id)
{
	return mesh->elem_get_area(mesh->msh, id);
}

double nb_mesh2D_elem_get_radius(const nb_mesh2D_t *mesh, uint32_t id)
{
	return mesh->elem_get_radius(mesh->msh, id);
}

double nb_mesh2D_elem_get_apotem(const nb_mesh2D_t *mesh, uint32_t id)
{
	return mesh->elem_get_apotem(mesh->msh, id);
}

uint32_t nb_mesh2D_elem_find_edge(const nb_mesh2D_t *mesh, uint32_t id,
				     uint16_t local_face_id)
{
	return mesh->elem_find_edge(mesh->msh, id, local_face_id);
}

double nb_mesh2D_elem_face_get_length(const nb_mesh2D_t *mesh,
					 uint32_t elem_id, uint16_t face_id)
{
	return mesh->elem_face_get_length(mesh->msh, elem_id, face_id);
}

double nb_mesh2D_elem_face_get_normal(const nb_mesh2D_t *mesh,
					 uint32_t elem_id, uint16_t face_id,
					 double normal[2])
{
	return mesh->elem_face_get_normal(mesh->msh, elem_id,
					  face_id, normal);
}

uint32_t nb_mesh2D_elem_face_get_left_ngb(const nb_mesh2D_t *mesh,
					     uint32_t elem_id,
					     uint16_t face_id)
{
	uint16_t N_adj = nb_mesh2D_elem_get_N_adj(mesh, elem_id);

	uint32_t left_elem;
	if (face_id < N_adj) {
		face_id = (face_id + 1) % N_adj;
		left_elem = nb_mesh2D_elem_get_ngb(mesh, elem_id,
						      face_id);
	} else {
		left_elem = nb_mesh2D_get_N_elems(mesh);
	}
	return left_elem;	
}

uint32_t nb_mesh2D_elem_face_get_right_ngb(const nb_mesh2D_t *mesh,
					      uint32_t elem_id,
					      uint16_t face_id)
{
	uint16_t N_adj = nb_mesh2D_elem_get_N_adj(mesh, elem_id);

	uint32_t right_elem;
	if (face_id < N_adj) {
		face_id = (face_id == 0)?(N_adj-1):(face_id-1);
		right_elem = nb_mesh2D_elem_get_ngb(mesh, elem_id,
						       face_id);
	} else {
		right_elem = nb_mesh2D_get_N_elems(mesh);
	}
	return right_elem;	
}

double nb_mesh2D_elem_ngb_get_normal(const nb_mesh2D_t *mesh,
					uint32_t elem_id, uint16_t ngb_id,
					double normal[2])
{
	return mesh->elem_ngb_get_normal(mesh->msh, elem_id,
					 ngb_id, normal);
}

uint16_t nb_mesh2D_elem_ngb_get_face(const nb_mesh2D_t *mesh,
					uint32_t elem_id, uint32_t ngb_id)
{
	uint16_t N_adj = nb_mesh2D_elem_get_N_adj(mesh, elem_id);
	uint16_t face_id = N_adj;
	for (uint16_t i = 0; i < N_adj; i++) {
		uint32_t ingb = nb_mesh2D_elem_get_ngb(mesh, elem_id, i);
		if (ingb == ngb_id) {
			face_id = i;
			break;
		}
	}
	return face_id;
}

double nb_mesh2D_elem_get_x(const nb_mesh2D_t *mesh, uint32_t id)
{
	return mesh->elem_get_x(mesh->msh, id);
}

double nb_mesh2D_elem_get_y(const nb_mesh2D_t *mesh, uint32_t id)
{
	return mesh->elem_get_y(mesh->msh, id);
}

uint32_t nb_mesh2D_elem_get_N_adj(const nb_mesh2D_t *mesh, uint32_t id)
{
	return mesh->elem_get_N_adj(mesh->msh, id);
}

uint32_t nb_mesh2D_elem_get_adj(const nb_mesh2D_t *mesh,
				   uint32_t elem_id, uint8_t adj_id)
{
	return mesh->elem_get_adj(mesh->msh, elem_id, adj_id);
}

uint32_t nb_mesh2D_elem_get_ngb(const nb_mesh2D_t *mesh,
				   uint32_t elem_id, uint8_t ngb_id)
{
	return mesh->elem_get_ngb(mesh->msh, elem_id, ngb_id);
}

bool nb_mesh2D_elem_has_ngb(const nb_mesh2D_t *mesh, uint32_t elem_id,
			       uint16_t face_id)
{
	return mesh->elem_has_ngb(mesh->msh, elem_id, face_id);
}

bool nb_mesh2D_elem_is_boundary(const nb_mesh2D_t *mesh, uint32_t elem_id)
{
	return mesh->elem_is_boundary(mesh->msh, elem_id);
}

uint32_t nb_mesh2D_get_invtx(const nb_mesh2D_t *mesh, uint32_t id)
{
	return mesh->get_invtx(mesh->msh, id);
}

uint32_t nb_mesh2D_insgm_get_N_nodes(const nb_mesh2D_t *mesh,
					  uint32_t id)
{
	return mesh->insgm_get_N_nodes(mesh->msh, id);
}

uint32_t nb_mesh2D_insgm_get_N_subsgm(const nb_mesh2D_t *mesh,
					 uint32_t id)
{
	uint32_t N = mesh->insgm_get_N_nodes(mesh->msh, id);
	return N - 1;
}

uint32_t nb_mesh2D_insgm_get_node(const nb_mesh2D_t *mesh,
				       uint32_t sgm_id, uint32_t node_id)
{
	return mesh->insgm_get_node(mesh->msh, sgm_id, node_id);
}

double nb_mesh2D_insgm_get_length(const nb_mesh2D_t *mesh,
				     uint32_t sgm_id)
{
	uint32_t last_vtx = mesh->insgm_get_N_nodes(mesh->msh, sgm_id) - 1;
	uint32_t n1 = mesh->insgm_get_node(mesh->msh, sgm_id, 0);
	uint32_t n2 = mesh->insgm_get_node(mesh->msh, sgm_id, last_vtx);
	double x1 = mesh->node_get_x(mesh->msh, n1);
	double y1 = mesh->node_get_y(mesh->msh, n1);
	double x2 = mesh->node_get_x(mesh->msh, n2);
	double y2 = mesh->node_get_y(mesh->msh, n2);
	return sqrt(POW2(x1 - x2) + POW2(y1 - y2));
}

double nb_mesh2D_insgm_subsgm_get_length(const nb_mesh2D_t *mesh,
					    uint32_t sgm_id,
					    uint32_t subsgm_id)
{
	uint32_t n1 = mesh->insgm_get_node(mesh->msh, sgm_id, subsgm_id);
	uint32_t n2 = mesh->insgm_get_node(mesh->msh, sgm_id, subsgm_id + 1);
	double x1 = mesh->node_get_x(mesh->msh, n1);
	double y1 = mesh->node_get_y(mesh->msh, n1);
	double x2 = mesh->node_get_x(mesh->msh, n2);
	double y2 = mesh->node_get_y(mesh->msh, n2);
	return sqrt(POW2(x1 - x2) + POW2(y1 - y2));
}

void nb_mesh2D_load_from_tessellator2D(nb_mesh2D_t *mesh,
			      nb_tessellator2D_t *t2D)
{
	mesh->load_from_mesh(mesh->msh, t2D);
}

void nb_mesh2D_set_nodal_permutation(nb_mesh2D_t *mesh,
					const uint32_t *perm)
{
	mesh->set_nodal_permutation(mesh->msh, perm);
}

void nb_mesh2D_get_enveloping_box(const nb_mesh2D_t *mesh,
				     double box[4])
{
	mesh->get_enveloping_box(mesh->msh, box);
}

bool nb_mesh2D_is_vtx_inside(const nb_mesh2D_t *mesh,
				double x, double y)
{
	return mesh->is_vtx_inside(mesh->msh, x, y);
}

void nb_mesh2D_build_model(const nb_mesh2D_t *mesh, nb_model_t *model)
{
	mesh->build_model(mesh->msh, model);
}

void nb_mesh2D_build_model_disabled_elems
			(const nb_mesh2D_t *mesh,
			 const bool *elems_enabled,
			 nb_model_t *model,
			 uint32_t *N_input_vtx,
			 uint32_t **input_vtx)
{
	mesh->build_model_disabled_elems(mesh->msh,
					 elems_enabled,
					 model,
					 N_input_vtx,
					 input_vtx);
}
