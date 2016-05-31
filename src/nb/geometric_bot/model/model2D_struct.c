#include <stdlib.h>

#include "nb/geometric_bot/model/model2D_struct.h"
#include "nb/geometric_bot/model/model2D.h"

inline void nb_model_alloc_vertices(nb_model_t *model)
{
	if (0 < model->N)
		model->vertex = calloc(2 * model->N, sizeof(*(model->vertex)));
}

inline void nb_model_alloc_edges(nb_model_t *model)
{
	if (0 < model->M)
		model->edge = calloc(2 * model->M, sizeof(*(model->edge)));
}

inline void nb_model_alloc_holes(nb_model_t *model)
{
	if (0 < model->H)
		model->holes = calloc(2 * model->H, sizeof(*(model->holes)));
}
