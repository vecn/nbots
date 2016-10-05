#ifndef __NB_GEOMETRIC_BOT_MESH_PARTITION_STRUCT_H__
#define __NB_GEOMETRIC_BOT_MESH_PARTITION_STRUCT_H__

#include <stdint.h>
#include <stdbool.h>

#include "nb/graphics_bot.h"

#include "nb/geometric_bot/model/model3D.h"

typedef struct {
	void (*draw_wires)(const void *msh,
			   nb_graphics_context_t *g);
	void (*draw_boundaries)(const void *msh,
				nb_graphics_context_t *g);
	void (*fill_elems)(const void *msh,
			   nb_graphics_context_t *g);
	void (*fill_elems_field_on_nodes)(const void *msh,
					  nb_graphics_context_t *g,
					  const double *normalized_field,
					  nb_graphics_palette_preset palette);
	void (*fill_elems_field_on_elems)(const void *msh,
					  nb_graphics_context_t *g,
					  const double *normalized_field,
					  nb_graphics_palette_preset palette);
	void (*fill_elems_classes)(const void *msh,
				   nb_graphics_context_t *g,
				   const uint8_t *class, uint8_t N_colors,
				   const nb_graphics_color_t *colors);
	void (*fill_nodes)(const void *msh,
			   nb_graphics_context_t *g);
	void (*fill_nodes_classes)(const void *msh,
				   nb_graphics_context_t *g,
				   const uint8_t *class, uint8_t N_colors,
				   const nb_graphics_color_t *colors);
} nb_partition_graphics_i;

struct nb_partition_s {
	void *msh;
	nb_partition_type type;

	void (*init)(void *msh);
	void (*finish)(void *msh);
	void (*copy)(void *msh, const void *mshsrc);
	void (*clear)(void *msh);

	uint32_t (*get_N_invtx)(const void *msh);
	uint32_t (*get_N_insgm)(const void *msh);
	uint32_t (*get_N_nodes)(const void *msh);
	uint32_t (*get_N_edges)(const void *msh);
	uint32_t (*get_N_elems)(const void *msh);
	double (*node_get_x)(const void *msh, uint32_t id);
	double (*node_get_y)(const void *msh, uint32_t id);
	uint32_t (*edge_get_1n)(const void *msh, uint32_t id);
	uint32_t (*edge_get_2n)(const void *msh, uint32_t id);
	void (*edge_get_midpoint)(const void *msh,
				  uint32_t face_id, double w,
				  double midpoint[2]);
	double (*edge_get_normal)(const void *msh, uint32_t face_id,
				  double normal[2]);
	double (*elem_get_x)(const void *msh, uint32_t id);
	double (*elem_get_y)(const void *msh, uint32_t id);
	double (*elem_get_area)(const void *msh, uint32_t id);
	double (*elem_get_radius)(const void *msh, uint32_t id);
	double (*elem_get_apotem)(const void *msh, uint32_t id);
	uint32_t (*elem_find_edge)(const void *msh, uint32_t id,
				   uint16_t local_face_id);
	double (*elem_face_get_length)(const void *msh,
				       uint32_t elem_id,
				       uint16_t face_id);
	double (*elem_face_get_normal)(const void *msh,
				       uint32_t elem_id, uint16_t face_id,
				       double normal[2]);
	double (*elem_ngb_get_normal)(const void *msh,
				      uint32_t elem_id, uint16_t ngb_id,
				      double normal[2]);
	uint32_t (*elem_get_N_adj)(const void *msh, uint32_t id);
	uint32_t (*elem_get_adj)(const void *msh,
				 uint32_t elem_id, uint8_t adj_id);
	uint32_t (*elem_get_ngb)(const void *msh,
				 uint32_t elem_id, uint8_t ngb_id);
	bool (*elem_has_ngb)(const void *msh, uint32_t elem_id,
			     uint16_t face_id);
	uint32_t (*get_invtx)(const void *msh, uint32_t id);
	uint32_t (*insgm_get_N_nodes)(const void *msh, uint32_t id);
	uint32_t (*insgm_get_node)(const void *msh, uint32_t sgm_id,
				     uint32_t node_id);
	void (*load_from_mesh)(void *msh, nb_mesh_t *mesh);
	void (*set_nodal_permutation)(void *msh, const uint32_t *perm);
	void (*get_enveloping_box)(const void *msh, double box[4]);
	bool (*is_vtx_inside)(const void *msh, double x, double y);
	double (*distort_with_field)(void *msh, 
				     nb_partition_entity field_entity,
				     double *disp,
				     double max_disp);
	void (*build_model)(const void *msh, nb_model_t *model);
	void (*build_model_disabled_elems)(const void *msh,
					   const bool *elems_enabled,
					   nb_model_t *model,
					   uint32_t *N_input_vtx,
					   uint32_t **input_vtx);
	nb_partition_graphics_i graphics;
};

#endif
