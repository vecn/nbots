#ifndef __NB_GEOMETRIC_BOT_MODEL_MODEL3D_H__
#define __NB_GEOMETRIC_BOT_MODEL_MODEL3D_H__

#include <stdbool.h>
#include <stdint.h>

#include "nb/graph_bot.h"

#include "nb/geometric_bot/model/model3D_struct.h"

uint16_t nb_model3D_get_memsize(void);
void nb_model3D_init(void *model_ptr);
void nb_model3D_copy(void *model_ptr, const void *src_model_ptr);
void nb_model3D_finish(void *model_ptr);
void nb_model3D_clear(void *model_ptr);

void nb_model3D_load(void *model_ptr, const char* filename);
void nb_model3D_load_box(void *model_ptr,
			 double x_min, double y_min, double z_min,
			 double x_max, double y_max, double z_max);
int nb_model3D_save(const void *model, const char* filename);
void nb_model3D_get_enveloping_box(const void *const model,
				   double box[6]);
uint32_t nb_model3D_get_N_vtx(const void *const model);
uint32_t nb_model3D_get_N_edges(const void *const model);
uint32_t nb_model3D_get_N_faces(const void *const model);

#endif
