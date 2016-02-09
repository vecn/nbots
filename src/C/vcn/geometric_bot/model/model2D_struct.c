#include <stdlib.h>

#include "model2D_struct.h"
#include "vcn/geometric_bot/model/model2D.h"

inline void model_alloc_vertices(vcn_model_t *model)
{
	if (0 < model->N)
		model->vertex = calloc(2 * model->N, sizeof(*(model->vertex)));
}

inline void model_alloc_edges(vcn_model_t *model)
{
	if (0 < model->M)
		model->edge = calloc(2 * model->M, sizeof(*(model->edge)));
}

inline void model_alloc_holes(vcn_model_t *model)
{
	if (0 < model->H)
		model->holes = calloc(2 * model->H, sizeof(*(model->holes)));
}
