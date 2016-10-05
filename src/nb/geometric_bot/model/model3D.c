/* Por implementar, responsable: Jorge Lopez */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>

#include "nb/geometric_bot/model/model3D_struct.h"
#include "nb/geometric_bot/model/model3D.h"

uint16_t nb_model3D_get_memsize(void)
{
	return 0;
}

void nb_model3D_init(void *model_ptr)
{
	;
}

void nb_model3D_copy(void *model_ptr, const void *src_model_ptr)
{
	;
}

void nb_model3D_finish(void *model_ptr)
{
	;
}

void nb_model3D_clear(void *model_ptr)
{
	;
}

void nb_model3D_load(void *model_ptr, const char* filename)
{
	;
}

void nb_model3D_load_box(void *model_ptr,
			 double x_min, double y_min, double z_min,
			 double x_max, double y_max, double z_max)
{
	;
}

int nb_model3D_save(const void *model, const char* filename)
{
	return 1;
}

void nb_model3D_get_enveloping_box(const void *const model, double box[6])
{
	;
}

uint32_t nb_model3D_get_N_vtx(const void *const model)
{
	;
}

uint32_t nb_model3D_get_N_edges(const void *const model)
{
	;
}

uint32_t nb_model3D_get_N_faces(const void *const model)
{
	;
}
