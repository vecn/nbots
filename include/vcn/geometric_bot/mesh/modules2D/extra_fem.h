#ifndef __VCN_GEOMETRIC_BOT_MESH_MODULES2D_EXTRA_FEM_H__
#define __VCN_GEOMETRIC_BOT_MESH_MODULES2D_EXTRA_FEM_H__

#ifdef __cplusplus
extern "C" {
#endif

	/**
	 * @brief Duplicate vertices which are the only intersection of two
	 * disjoint portions of the mesh in order to disconnect such portions.
	 *
	 * @param[in] mesh Mesh to be processed.
	 */
	void vcn_mesh_duplicate_one_point_connections(vcn_mesh_t* mesh);


#ifdef __cplusplus
}
#endif

#endif
