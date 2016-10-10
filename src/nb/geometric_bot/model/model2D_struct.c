#include <stdlib.h>

#include "nb/memory_bot.h"
#include "nb/geometric_bot/model/model2D_struct.h"
#include "nb/geometric_bot/model/model2D.h"

inline void nb_model_alloc_vertices(nb_model_t *model)
{
	if (0 < model->N) {
		model->vertex = nb_allocate_zero_mem(2 * model->N * 
						     sizeof(*(model->vertex)));
	}
}

inline void nb_model_alloc_edges(nb_model_t *model)
{
	if (0 < model->M) {
		model->edge = nb_allocate_zero_mem(2 * model->M *
						   sizeof(*(model->edge)));
	}
}

inline void nb_model_alloc_holes(nb_model_t *model)
{
	if (0 < model->H) {
		model->holes = nb_allocate_zero_mem(2 * model->H *
						    sizeof(*(model->holes)));
	}
}
