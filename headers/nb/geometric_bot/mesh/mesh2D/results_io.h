#ifndef __NB_GEOMETRIC_BOT_MESH_MESH2D_RESULTS_IO_H__
#define __NB_GEOMETRIC_BOT_MESH_MESH2D_RESULTS_IO_H__


int nb_mesh2D_read_and_draw_results(const char *dir_saved_results,
				    const char *dir_output,
				    /* Show progress can be NULL */
				    void (*show_progress)(float prog));

#endif
