#ifndef __NB_PDE_BOT_CONTROL_VOLUME_INTEGRATION_MESH_H__
#define __NB_PDE_BOT_CONTROL_VOLUME_INTEGRATION_MESH_H__

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#include "nb/solver_bot.h"
#include "nb/geometric_bot.h"

typedef struct subface_s subface_t;

typedef struct {
	double nf[2];
	uint32_t elems[2];
	double x1[2], x2[2];
	uint16_t N_sf;
	subface_t **subfaces;
} face_t;

struct subface_s {
	uint8_t N_int;/*     Zero: Pairwise      */
                      /* Not zero: Simplex-wise  */
	double x1[2], x2[2];
	uint32_t trg_id;
};

uint32_t nb_cvfa_get_integration_mesh_memsize(void);
void nb_cvfa_init_integration_mesh(nb_mesh2D_t *intmsh);
void nb_cvfa_load_integration_mesh(nb_mesh2D_t *intmsh, uint32_t N,
				   const double *xc);

void nb_cvfa_correlate_mesh_and_integration_mesh
					(const nb_mesh2D_t *part,
					 const nb_mesh2D_t *intmsh,
					 nb_graph_t *trg_x_vol);

void nb_cvfa_get_adj_graph(const nb_mesh2D_t *intmsh,
			   const nb_graph_t *trg_x_vol,
			   nb_graph_t *graph);

void nb_cvfa_load_trg_points(const nb_mesh2D_t *intmsh,
			     uint32_t trg_id, double t1[2],
			     double t2[2], double t3[2]);
void nb_cvfa_init_global_matrix(nb_sparse_t **K, const nb_graph_t *trg_x_vol,
				const nb_mesh2D_t *intmsh);

void nb_cvfa_load_faces(const nb_mesh2D_t *mesh,
			const nb_mesh2D_t *intmsh,
			const nb_graph_t *trg_x_vol,
			face_t **faces);
void nb_cvfa_finish_faces(uint32_t N_faces, face_t **faces);

#endif
