#ifndef __NB_GEOMETRIC_BOT_MESH_PARTITION_PRIVATE_H__
#define __NB_GEOMETRIC_BOT_MESH_PARTITION_PRIVATE_H__

#include <stdint.h>

#include "partition_struct.h"

typedef struct {
	void (*node_move_x)(void *msh, uint32_t i, double u);
	void (*node_move_y)(void *msh, uint32_t i, double v);
	void (*elem_move_x)(void *msh, uint32_t i, double u);
	void (*elem_move_y)(void *msh, uint32_t i, double v);
} nb_partition_private_i;

void nb_partition_init_private_interface(nb_partition_private_i *priv,
					 const nb_partition_t *part);

#endif
