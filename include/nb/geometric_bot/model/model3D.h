#ifndef __NB_GEOMETRIC_BOT_MODEL_MODEL3D_H__
#define __NB_GEOMETRIC_BOT_MODEL_MODEL3D_H__

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <time.h>

#include "nb/graph_bot.h"

#include "nb/geometric_bot/model/model3D_struct.h"

uint64_t nb_model3D_get_memsize(uint32_t n_triangles);
void nb_model3D_init(void *model_ptr);
void nb_model3D_copy(void *model_ptr, const void *src_model_ptr);
void nb_model3D_finish(void *model_ptr);
void nb_model3D_clear(void *model_ptr);

void nb_model3D_load(void *model_ptr, const char* filename, uint32_t n_triangles);
void nb_model3D_load_box(void *model_ptr,
			 double x_min, double y_min, double z_min,
			 double x_max, double y_max, double z_max);
int nb_model3D_save(const void *model, const char* filename);
void nb_model3D_get_enveloping_box(const void *const model,
				   double box[6]);
uint32_t nb_model3D_get_N_vtx(const void *const model);
uint32_t nb_model3D_get_N_edges(const void *const model);
uint32_t nb_model3D_get_N_face(const void *const model);


uint32_t nb_model3D_get_number_of_triangles_in_ASCII_STL( const char *name );
void nb_model3D_read_header_ASCII_STL( FILE *fp , char *header );
nb_model3D_t* nb_model3D_assign_memory_on_struct_for_ASCII_STL(uint64_t memsize , 
                                                               uint32_t n_triangles);
void nb_model3D_load_model_from_ASCII_STL( const char *name , nb_model3D_t *boundary_t ,
                                            uint32_t n_triangles);
void nb_model3D_set_nf(int position,double value, void *boundary_t);
void nb_model3D_set_vtx(int position,double value, void *boundary_t);
void nb_model3D_set_adj(int position,uint32_t value, void *boundary_t);
void nb_model3D_set_N_vtx(uint32_t value, void *boundary_t);
void nb_model3D_set_N_face(uint32_t value, void *boundary_t);
double nb_model3D_get_nf(int position, void *boundary_t);
double nb_model3D_get_vtx(int position, void *boundary_t);
uint32_t nb_model3D_get_adj(int position, void *boundary_t);
#endif
