#include <stdbool.h>
#include <stdint.h>

#include "nb/memory_bot/nb_malloc.h"
#include "nb/memory_bot/membank.h"

#define NB_MEMBANK_STATIC_SIZE 5000
#define NB_MEMBANK_MASK_SIZE 625 /* 5000 / 8 */

struct nb_membank_s {
	uint16_t type_size;

	uint16_t N_max_static;
	uint16_t N_static;
	uint16_t N_static_free;
	char static_mem[NB_MEMBANK_STATIC_SIZE];
	char static_mask[NB_MEMBANK_MASK_SIZE];
	
	uint16_t N_max_dynamic;
	uint16_t N_dynamic;
	uint16_t N_dynamic_free;
	char *dynamic_mask;
	char *dynamic_mask;
};

static void set_mask(char mask[NB_MEMBANK_MASK_SIZE], int i);
static bool get_mask(char mask[NB_MEMBANK_MASK_SIZE], int i);
static void unset_mask(char mask[NB_MEMBANK_MASK_SIZE], int i);

void nb_membank_init(nb_membank_t *membank, uint16_t type_size)
{
	memset(membank, 0, sizeof(*membank));
	membank->type_size = type_size;
	membank->N_max_static = NB_MEMBANK_STATIC_SIZE / type_size;
}

void* nb_membank_calloc(nb_membank_t *membank)
{
	void *mem;
	if (membank->N_static < membank->N_max_static) {
		mem = membank->static_mem +
			membank->N_static * membank->type_size;
		set_mask(membank->static_mask, membank->N_static);
		membank->N_static += 1;
	} else if (membank->N_static_free > membank->N_max_static / 3) {
		int i = 0;
		while (get_mask(membank->static_mask, i))
			i += 1;
		mem = membank->static_mem + i * membank->type_size;
		set_mask(membank->static_mask, i);
		membank->N_static_free -= 1;
	} else {
		/* AQUI VOY */
		mem = nb_calloc(membank->type_size);
	}
	return mem;
}

static void set_mask(char mask[NB_MEMBANK_MASK_SIZE], int i)
{

}

static bool get_mask(char mask[NB_MEMBANK_MASK_SIZE], int i)
{

}

void nb_membank_free(nb_membank_t *membank, void *obj)
{
	char *mem = obj;
	if (mem >= membank->static_mem && 
	    mem < membank->static_mem + NB_MEMBANK_STATIC_SIZE) {
		int i = (mem - membank->static_mem) / membank->type_size;
		unset_mask(membank->static_mask, i);
		membank->N_static_free += 1;
		if (i == membank->N_static) {
			do {
				membank->N_static -= 1;
			} while (!get_mask(membank->static_mask,
					   membank->N_static));
		}
	} else {
		/* AQUI VOY */
		nb_free(mem);
	}
}

static void unset_mask(char mask[NB_MEMBANK_MASK_SIZE], int i)
{

}
