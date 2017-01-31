#include <stdlib.h>
#include <string.h>

#include "nb/memory_bot.h"
#include "nb/geometric_bot/model/model2D_struct.h"
#include "nb/geometric_bot/model/model2D.h"

void nb_model_alloc_vertices(nb_model_t *model)
{
	if (0 < model->N) {
		model->vertex = nb_allocate_zero_mem(2 * model->N * 
						     sizeof(*(model->vertex)));
	}
}

void nb_model_alloc_edges(nb_model_t *model)
{
	if (0 < model->M) {
		model->edge = nb_allocate_zero_mem(2 * model->M *
						   sizeof(*(model->edge)));
	}
}

void nb_model_alloc_holes(nb_model_t *model)
{
	if (0 < model->H) {
		model->holes = nb_allocate_zero_mem(2 * model->H *
						    sizeof(*(model->holes)));
	}
}

void nb_model_load_from_arrays(nb_model_t *model,
			       uint32_t N_vtx, double *vtx,
			       uint32_t N_sgm, uint32_t *sgm,
			       uint32_t N_holes, double *holes)
{
	nb_model_clear(model);
	model->N = N_vtx;
	nb_model_alloc_vertices(model);
	memcpy(model->vertex, vtx, 2 * N_vtx * sizeof(*(model->vertex)));
	
	model->M = N_sgm;
	if (N_sgm > 0) {
		nb_model_alloc_edges(model);
		memcpy(model->edge, sgm,
		       2 * N_sgm * sizeof(*(model->edge)));
	}

	model->H = N_holes;
	if (N_holes > 0) {
		nb_model_alloc_holes(model);
		memcpy(model->holes, holes,
		       2 * N_holes * sizeof(*(model->holes)));
	}
}
