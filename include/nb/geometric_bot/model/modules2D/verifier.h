#ifndef __NB_GEOMETRIC_BOT_MODEL_MODULES2D_VERIFIER_H__
#define __NB_GEOMETRIC_BOT_MODEL_MODULES2D_VERIFIER_H__

#include <stdbool.h>
#include <stdint.h>

#include "nb/geometric_bot/model/model2D.h"
  
/**
 * @brief
 * @param[in] model Model to be verified.
 * @param[out] ids_causing_error NULL if not required. Stores the ids of 
 * the model elements producing the inconsistency, such elements could be
 * vertices or segments.
 * @return Error code:
 * + <b>0</b> No error. The model could be meshed.
 * + <b>1</b> There are no vertices.
 * + <b>2</b> There are no segments.
 * + <b>3</b> There are repeated vertices.
 * + <b>4</b> Vertex id (forming an edge) out of bounds.
 * + <b>5</b> There are repeated segments.
 * + <b>6</b> There are intersecting segments.
 * + <b>7</b> There is a vertex upon a segment.
 * + <b>8</b> Unknown error (Unlikely but possible due to numerical error).
 * + <b>9</b> Unclosed shape.
 */
int vcn_model_verify_consistence(const vcn_model_t *const model,
				 uint32_t ids_causing_error[2]);

bool vcn_model_have_vertices(const vcn_model_t *const model);
bool vcn_model_have_edges(const vcn_model_t *const model);
bool vcn_model_have_repeated_vertices(const vcn_model_t *const model,
				      uint32_t repeated_ids[2]);
bool vcn_model_have_incoherent_edges(const vcn_model_t *const model,
				     uint32_t ids_edge_and_vtx[2]);
bool vcn_model_have_repeated_edges(const vcn_model_t *const model,
				   uint32_t repeated_ids[2]);
bool vcn_model_have_intersected_edges(const vcn_model_t *const model,
				      uint32_t intersected_ids[2]);
bool vcn_model_have_vtx_intersecting_edges(const vcn_model_t *const model,
					   uint32_t ids_edge_and_vtx[2]);
bool vcn_model_have_unclosed_boundary(const vcn_model_t *const model);
bool vcn_model_is_continuum(const vcn_model_t *model);

uint16_t vcn_model_get_N_subareas(const vcn_model_t *model);


#endif
