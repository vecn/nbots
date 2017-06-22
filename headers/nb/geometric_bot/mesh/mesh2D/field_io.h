#ifndef __NB_GEOMETRIC_BOT_MESH_MESH2D_FIELD_IO_H__
#define __NB_GEOMETRIC_BOT_MESH_MESH2D_FIELD_IO_H__

#include "nb/geometric_bot.h"

int nb_mesh2D_field_read_and_draw(const char *dir_saved_results,
				  const char *dir_output,
				  /* Show progress can be NULL */
				  void (*show_progress)(float prog));

int nb_mesh2D_field_get_N_steps(const char *filename, uint32_t *N_steps);

int nb_mesh2D_field_read_step_field(const char *filename, uint32_t step,
				    double *time, const char *field_name,
				    double *field);

#endif
