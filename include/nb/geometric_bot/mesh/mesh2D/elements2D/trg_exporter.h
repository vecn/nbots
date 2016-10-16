#ifndef __NB_GEOMETRIC_BOT_MESH_MESH2D_ELEMENTS2D_TRG_EXPORTER_H__
#define __NB_GEOMETRIC_BOT_MESH_MESH2D_ELEMENTS2D_TRG_EXPORTER_H__

#include <stdbool.h>
#include <stdint.h>
#include "nb/geometric_bot/mesh/tessellator2D.h"
#include "nb/graph_bot.h"

typedef struct {
	void *structure;
	/* Mandatory */
	void (*set_N_vtx)(void*, uint32_t N);
	void (*allocate_vtx)(void*);
	void (*start_vtx_access)(void*);
	void (*set_vtx)(void*, uint32_t i, double x, double y);
	void (*stop_vtx_access)(void*);

	/* Optional (can be NULL) */
	void (*set_N_edg)(void*, uint32_t N); /* Used as flag */
	void (*allocate_edg)(void*);
	void (*start_edg_access)(void*);
	void (*set_edg)(void*, uint32_t i, uint32_t v1, uint32_t v2);
	void (*stop_edg_access)(void*);

	/* Optional (can be NULL) */
	void (*set_N_trg)(void*, uint32_t N); /* Used as flag */
	void (*allocate_trg)(void*, bool include_neighbours);
	void (*start_trg_access)(void*);
	void (*set_trg)(void*, uint32_t i, uint32_t v1,
			uint32_t v2, uint32_t v3);
	void (*stop_trg_access)(void*);

	/* Optional (can be NULL) */
	void (*start_trg_neighbours_access)(void*);
	void (*set_trg_neighbours)(void*, uint32_t i, uint32_t t1,
				   uint32_t t2, uint32_t t3); /* Used as flag */
	void (*stop_trg_neighbours_access)(void*);

	/* Optional (can be NULL) */
	void (*set_N_input_vtx)(void*, uint32_t N);/* Used as flag */
	void (*allocate_input_vtx)(void*);
	void (*start_input_vtx_access)(void*);
	void (*set_input_vtx)(void*, uint32_t i, uint32_t vtx_id);
	void (*stop_input_vtx_access)(void*);

	/* Optional (can be NULL) */
	void (*set_N_input_sgm)(void*, uint32_t N);/* Used as flag */
	void (*allocate_input_sgm_table)(void*);
	void (*start_input_sgm_table_access)(void*);
	void (*input_sgm_set_N_vtx)(void*, uint32_t i, uint32_t N);
	void (*input_sgm_allocate_vtx)(void*, uint32_t i);
	void (*input_sgm_start_access)(void*, uint32_t isgm);
	void (*input_sgm_set_vtx)(void*, uint32_t ivtx, uint32_t vtx_id);
	void (*input_sgm_stop_access)(void*);
	void (*stop_input_sgm_table_access)(void*);
} nb_trg_exporter_interface_t;

void nb_tessellator2D_export(const nb_tessellator2D_t *const mesh,
		     nb_trg_exporter_interface_t *exp);


#endif
