#ifndef __NB_MEMORY_BOT_MEMBANK_H__
#define __NB_MEMORY_BOT_MEMBANK_H__

#include <stdint.h>

typedef struct nb_membank_s nb_membank_t;

uint32_t nb_membank_get_memsize(void);
void nb_membank_init(nb_membank_t *membank, uint16_t type_size);
void* nb_membank_allocate_mem(nb_membank_t *membank);
void nb_membank_free_mem(nb_membank_t *membank, void *obj);
void nb_membank_finish(nb_membank_t *membank);
void nb_membank_clear(nb_membank_t *membank);
void nb_membank_merge(nb_membank_t *membank, nb_membank_t *extension);
void nb_membank_set_N_x_block(nb_membank_t *membank, uint16_t N_x_block);

#endif
