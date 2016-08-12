#ifndef __NB_GEOMETRIC_BOT_MESH2D_DISKS_STRUCT_H__
#define __NB_GEOMETRIC_BOT_MESH2D_DISKS_STRUCT_H__

/**
 * @brief Read-only mesh structure, which stores a 2D sphere packing
 * produced from the triangulation in vcn_mesh_t.
 */
typedef struct vcn_mshpack_s {
	uint32_t N_spheres;
	double* centers;
	double* radii;
	uint32_t* N_adj;
	uint32_t** adj;
} vcn_mshpack_t;

#endif
