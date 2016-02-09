#ifndef __VCN_GEOMETRIC_BOT_MESH_MESH2D_STRUCTURES_H__
#define __VCN_GEOMETRIC_BOT_MESH_MESH2D_STRUCTURES_H__

#include "vcn/geometric_bot/point2D.h"
#include "vcn/geometric_bot/knn/bins2D.h"

#define _VCN_INPUT_VTX ((void*)0x1)

typedef vcn_point2D_t msh_vtx_t;
typedef struct msh_edge_s msh_edge_t;
typedef struct msh_trg_s msh_trg_t;

typedef struct {
	uint32_t id;
	void* data;
} attr_t;

typedef struct {
	uint32_t idx;
	void* attr;
	msh_edge_t* prev;
	msh_edge_t* next;
} input_sgm_attr_t;

struct msh_edge_s {
	void* attr;
	msh_vtx_t *v1, *v2;
	msh_trg_t *t1, *t2;
};

struct msh_trg_s {
	/* Mesh Triangle 
	 *                    /\
	 *                   /  \
	 *                  /    \ 
	 *                 /__T3__\  
	 *                  ^ ^ |
	 *                  | | |
	 *                 _|_v_|_ s3 
	 *                  | ^ |
	 *                  | | |
	 *             V1 ._|_v_v_. V3
	 *    /\-----\---->\     /-----/---->/\
	 *   /  \-> <-\-> <-\   /-> <-/-> <-/  \
	 *  /_T1_\<----\-----\ /<----/-----/_T2_\ 
	 *             S1     °    S2 
	 *                    V2
	 *
	 */
	void *attr;
	msh_trg_t *t1, *t2, *t3;   /* Neighbouring triangles */
	msh_edge_t *s1, *s2, *s3;  /* Subsegments */
	msh_vtx_t *v1, *v2, *v3;   /* Vertex */
};

struct vcn_mesh_s {
	/* Arrays to relate the mesh with the PSLG input */
	uint32_t N_input_vtx;
	msh_vtx_t** input_vtx;
	uint32_t N_input_sgm;      
	msh_edge_t** input_sgm; /* Each input_sgm is the root of a double 
				 * linked list conatining the subsegments
				 * forming this input segment. */
	/* Data structures to handle mesh */
	vcn_bins2D_t* ug_vtx;       /* Grid to sort vertices */
	vcn_container_t* ht_trg;   /* Hash table of triangles */
	vcn_container_t* ht_edge;   /* Hash table of segments */

	/* Scale and shift values to handle floating point error */
	double scale;
	double xdisp;
	double ydisp;

	/* Size constraints */
	uint32_t max_vtx;
	uint32_t max_trg;

	/* Geometric constrains */
	double min_angle;
	double cr2se_ratio;
	double max_edge_length;
	double max_subsgm_length;

	/* Density data */
	double (*density)(const double x[2], const void *data);
	const void *density_data;

	/* Refinement algorithm type */
	int refiner_type;

	/* Tasks functions */
	void (*do_after_insert_trg)(const vcn_mesh_t *const);
	void (*do_after_insert_vtx)(const vcn_mesh_t *const);
};

bool medge_is_subsgm(const msh_edge_t *const sgm);
void medge_set_as_subsgm(msh_edge_t *const sgm,
			 uint32_t idx,
			 const msh_edge_t *const prev,
			 const msh_edge_t *const next);
void medge_update_subsgm_next(msh_edge_t *const sgm, 
			      const msh_edge_t *const next);
void medge_update_subsgm_prev(msh_edge_t *const sgm,
			      const msh_edge_t *const prev);
uint32_t medge_subsgm_idx(const msh_edge_t *const sgm);
msh_edge_t* medge_subsgm_prev(const msh_edge_t *const sgm);
msh_edge_t* medge_subsgm_next(const msh_edge_t *const sgm);
void medge_subsgm_set_attribute(msh_edge_t* sgm, void* attr);
void* medge_subsgm_get_attribute(const msh_edge_t* const sgm);
void medge_destroy_subsgm_attribute(msh_edge_t *const sgm);
void medge_set_length(msh_edge_t *const sgm);
void medge_update_length(msh_edge_t *const sgm);
double medge_get_computed_length(const msh_edge_t *const sgm);
void medge_destroy_length_attribute(msh_edge_t *const sgm);
bool medge_has_a_triangle_to_the_left(const msh_edge_t *const sgm, 
				      const msh_vtx_t *const v1);
void medge_connect_triangles(msh_edge_t *const restrict sgm);
void medge_connect_triangle(msh_edge_t *const sgm, 
			    msh_trg_t *const trg,
			    const msh_vtx_t *const v1, 
			    const msh_vtx_t *const v2);
msh_trg_t* medge_get_triangle_to_the_left
                           (const msh_edge_t *const sgm,
			    const msh_vtx_t *const v1);
msh_trg_t* medge_get_opposite_triangle(const msh_edge_t *const sgm,
				       const msh_trg_t *const trg);
msh_edge_t* medge_get_CW_subsgm(const msh_edge_t *const sgm,
				const msh_vtx_t *const vtx);
msh_edge_t* medge_get_CCW_subsgm(const msh_edge_t *const sgm,
				 const msh_vtx_t *const vtx);
void medge_flip_without_dealloc(msh_edge_t* shared_sgm);
msh_trg_t* mtrg_create(void);
void mtrg_set_quality_and_size(msh_trg_t *const trg);
double mtrg_get_computed_quality(const msh_trg_t *const trg);
double mtrg_get_computed_size(const msh_trg_t *const trg);
void mtrg_destroy_quality_and_size_attributes(msh_trg_t *const trg);
bool mtrg_has_an_input_vertex(const msh_trg_t *const trg);
msh_vtx_t* mtrg_get_opposite_vertex_guided(const msh_trg_t *const trg, 
					   const msh_edge_t *const sgm,
					   bool same_vertices_order);
msh_vtx_t* mtrg_get_opposite_vertex(const msh_trg_t *const trg, 
				    const msh_edge_t *const sgm);
msh_edge_t* mtrg_get_largest_edge(const msh_trg_t *const trg);
msh_edge_t* mtrg_get_shortest_edge(const msh_trg_t *const trg);
msh_edge_t* mtrg_get_opposite_edge(const msh_trg_t *const trg, 
				   const msh_vtx_t *const vtx);
msh_edge_t* mtrg_get_right_edge(const msh_trg_t *const trg, 
				const msh_vtx_t *const vtx);
msh_edge_t* mtrg_get_left_edge(const msh_trg_t *const trg, 
			       const msh_vtx_t *const vtx);
void mtrg_get_complement_edges(const msh_trg_t *const trg,
			       const msh_edge_t *const edge,
			       msh_edge_t* edge_complement[2]);
msh_trg_t* mtrg_get_right_triangle(const msh_trg_t *const trg, 
				   const msh_vtx_t *const vtx);
msh_trg_t* mtrg_get_left_triangle(const msh_trg_t *const trg, 
				  const msh_vtx_t *const vtx);
void mtrg_vanish_from_neighbour(const msh_trg_t *const trg, 
				msh_trg_t *const nb_trg);
void mtrg_disconnect(const msh_trg_t *const trg);
msh_edge_t* mtrg_get_CCW_edge(const msh_trg_t *const trg,
			      const msh_edge_t *const edge);
msh_edge_t* mtrg_get_CW_edge(const msh_trg_t *const trg,
			     const msh_edge_t *const edge);
msh_edge_t* mesh_insert_edge(vcn_container_t *const ht_edge,
			     const msh_vtx_t *const v1, 
			     const msh_vtx_t *const v2);
msh_edge_t* mesh_exist_edge_guided(vcn_container_t *const ht_edge,
				   const msh_vtx_t *const v1,
				   const msh_vtx_t *const v2);
msh_edge_t* mesh_exist_edge(vcn_container_t *const ht_edge,
			    const msh_vtx_t *const v1,
			    const msh_vtx_t *const v2);
void mesh_add_triangle(vcn_mesh_t *const mesh, msh_trg_t *const trg);
void mesh_substract_triangle(vcn_mesh_t *const mesh,
			     const msh_trg_t *const trg);
uint32_t hash_key_edge(const void *const  edge_ptr);
bool are_equal_edge(const void *const edge1_ptr, const void *const edge2_ptr);
msh_trg_t* mesh_locate_vtx(const vcn_mesh_t *const mesh,
			  const msh_vtx_t *const v);
void mesh_get_extern_scale_and_disp(const vcn_mesh_t *const mesh,
				    const double *const internal,
				    double external[2]);
#endif
