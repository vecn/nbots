#ifndef __NB_GEOMETRIC_BOT_MESH_ELEMENTS2D_TRIANGLES_STRUCT_H__
#define __NB_GEOMETRIC_BOT_MESH_ELEMENTS2D_TRIANGLES_STRUCT_H__

#define vcn_msh3trg_get_1st_vtx_id_of_trg(msh3trg, trg_id)\
  ((msh3trg)->vertices_forming_triangles[(trg_id) * 3])
#define vcn_msh3trg_get_2nd_vtx_id_of_trg(msh3trg, trg_id)\
  ((msh3trg)->vertices_forming_triangles[(trg_id)*3+1])
#define vcn_msh3trg_get_3rd_vtx_id_of_trg(msh3trg, trg_id)\
  ((msh3trg)->vertices_forming_triangles[(trg_id)*3+2])
#define vcn_msh3trg_view_vtx(msh3trg, vtx_id)\
  (&((msh3trg)->vertices[(vtx_id)*2]))

/**
 * @brief Read-only mesh structure, which stores the same mesh than
 * vcn_mesh_t (a Delaunay triangulation), but stored in easy-to-read
 * data structures.
 */
typedef struct vcn_msh3trg_s {
	uint32_t N_vertices;
	double *vertices;   /* Vertices coordinates concatenated */

	uint32_t N_edges;
	uint32_t *edges;

	uint32_t N_triangles;
	uint32_t *vertices_forming_triangles; /* Connectivity matrix */
	uint32_t *triangles_sharing_sides;    /* Neighbours */

	uint32_t N_input_vertices;
	uint32_t *input_vertices; /* Ids of vtx corresponding to input vtx */

	uint32_t N_input_segments;
	/* Number of subsegments forming the input segment */
	uint32_t *N_subsgm_x_inputsgm;
	/* Sequence of vertices ids forming the input segments ************/
	uint32_t** meshvtx_x_inputsgm;                                  /**/
	/* Important:                                                     */
	/*    The length of meshvtx_x_inputsgm[i] is equal to             */
	/*              N_subsgm_x_inputsgm[i] + 1,                       */
	/*    since each subsegment is characterized by the first         */
	/*    vertex, except the last in the chain, which needs           */
	/*    its both vertices to be defined.                            */
	/******************************************************************/
} vcn_msh3trg_t;


#endif