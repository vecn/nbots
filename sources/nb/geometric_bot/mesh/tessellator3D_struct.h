#ifndef __NB_GEOMETRIC_BOT_MESH_TESSELLATOR_3D_STRUCT_H__
#define __NB_GEOMETRIC_BOT_MESH_TESSELLATOR_3D_STRUCT_H__

#include <stdint.h>

struct nb_tessellator3D_s {
	void *t3D;
	nb_tessellator3D_type type;

	uint32_t (*get_memsize)(void);
	void (*init)(void *t3D);
	void (*copy)(void *t3D,  const void *src);
	void (*finish)(void *t3D);

	void (*clear)(void* t3D);
	void (*set_density)(void* t3D,
			    double (*density)(const double x[3],
					      const void *data),
			    const void *density_data);
	void (*unset_density)(void* t3D);
	bool (*is_empty)(const void *const t3D);
	void (*generate_from_model)(void *t3D,
				    const nb_model3D_t *const model);
	void (*get_simplest_from_model)(void *t3D,
					const nb_model3D_t *const  model);
	bool (*is_vtx_inside)(const void *const t3D,
			      const double vtx[3]);
	void (*refine)(void *t3D);
	bool (*insert_vtx)(void *t3D,
			   const double vtx[3]);
	uint32_t (*get_N_vtx)(const void *t3D);
	uint32_t (*get_N_edges)(const void *t3D);
	uint32_t (*get_N_faces)(const void *t3D);
	uint32_t (*get_N_tetra)(const void *t3D);
	double (*get_vol)(const void *t3D);
};

#endif
