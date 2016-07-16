#ifndef __NB_MEMORY_BOT_MEMBANK_H__
#define __NB_MEMORY_BOT_MEMBANK_H__

#include <stdint.h>

typedef struct nb_membank_s nb_membank_t;

uint32_t nb_membank_get_memsize(void);
void nb_membank_init(nb_membank_t *membank, uint16_t type_size);
void* nb_membank_calloc(nb_membank_t *membank);
void nb_membank_free(nb_membank_t *membank, void *obj);
void nb_membank_finish(nb_membank_t *membank);
void nb_membank_clear(nb_membank_t *membank);

#endif
