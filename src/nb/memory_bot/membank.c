#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

#include "nb/memory_bot/nb_malloc.h"
#include "nb/memory_bot/membank.h"

#define STATIC_SIZE 8000
#define MASK_SIZE 1000
#define DIV_FREED_TO_ALLOC 5
#define MIN_MAX_DYNAMIC 100

#define MAX(a,b) (((a)>(b))?(a):(b))

#define SET_MASK(mask, i)			\
	((mask)[(i)/8] |= 1 << ((i)%8))
#define IS_MASK_ENABLED(mask, i)		\
	((mask)[(i)/8] & (1 << ((i)%8)))
#define UNSET_MASK(mask, i)			\
	((mask)[(i)/8] &=  ~(1 << ((i)%8)))
#define IS_MASK_DISABLED(mask, i)		\
	(!IS_MASK_ENABLED(mask, i))
#define IS_MASK_BYTE_FULL(mask, char_id)	\
	(~((mask)[(char_id)]) == 0)
#define IS_MASK_BYTE_ENABLED(mask, char_id, bit_id)	\
	((mask)[(char_id)] & (1 << (bit_id)))

typedef struct block_s block_t;

struct block_s {
	uint16_t N_max;
	uint16_t N_align;
	uint16_t N_free;
	char *buffer;
	char *mask;
	block_t *next;
};

struct nb_membank_s {
	uint16_t type_size;
	uint16_t N_max_dynamic;
	block_t block;
	char static_buffer[STATIC_SIZE];
	char static_mask[MASK_SIZE];	
};

static void* block_calloc(block_t *block, uint16_t type_size);
static block_t *block_create(uint16_t type_size, uint16_t N_max);
static void set_backward_aligner(block_t *block, int i);
static bool block_is_in_buffer(const block_t *block,
				  uint16_t type_size, char *mem);
static void block_free(block_t *block, uint16_t type_size,
			  char *mem);

uint32_t nb_membank_get_memsize(void)
{
	return sizeof(nb_membank_t);
}

void nb_membank_init(nb_membank_t *membank, uint16_t type_size)
{
	memset(membank, 0, sizeof(*membank));
	membank->type_size = type_size;
	membank->N_max_dynamic = MAX(MIN_MAX_DYNAMIC,
				     STATIC_SIZE / type_size);
	membank->block.N_max = STATIC_SIZE / type_size;
	membank->block.N_free = STATIC_SIZE / type_size;
	membank->block.buffer = membank->static_buffer;
	membank->block.mask = membank->static_mask;
}

void *nb_membank_calloc(nb_membank_t *membank)
{
	block_t *block = &(membank->block);
	void *mem = NULL;
	while (NULL != block) {
		mem = block_calloc(block, membank->type_size);
		if (NULL != mem) {
			goto EXIT;
		} else {
			if (NULL == block->next)
				block->next = 
					block_create(membank->type_size,
						     membank->N_max_dynamic);
			block = block->next;
		}
	}
EXIT:
	return mem;
}

static void *block_calloc(block_t *block, uint16_t type_size)
{
	void *mem = NULL;
	int i;
	if (block->N_align < block->N_max) {
		i = block->N_align;
		block->N_align += 1;
	} else if (block->N_free > block->N_max / DIV_FREED_TO_ALLOC) {
		int char_id = 0;
		while (IS_MASK_BYTE_FULL(block->mask, char_id))
			char_id += 1;

		int bit_id = 0;
		while (IS_MASK_BYTE_ENABLED(block->mask, char_id, bit_id))
			bit_id += 1;

		i = char_id * 8 + bit_id;
	} else {
		goto EXIT;
	}
	mem = block->buffer + i * type_size;
	memset(mem, 0, type_size);
	SET_MASK(block->mask, i);
	block->N_free -= 1;
EXIT:
	return mem;
}

static block_t *block_create(uint16_t type_size, uint16_t N_max)
{
	uint32_t block_size = sizeof(block_t);
	uint32_t buffer_size = type_size * N_max;
	uint32_t mask_size = N_max / 8 + ((N_max % 8 > 0)?1:0);
	uint32_t size = block_size + buffer_size + mask_size;

	char *memory = nb_calloc(size);

	block_t *block = (void*) memory;
	block->N_max = N_max;
	block->N_free = N_max;
	block->buffer = memory + block_size;
	block->mask = memory + block_size + buffer_size;
	return block;
}


void nb_membank_free(nb_membank_t *membank, void *obj)
{
	block_t *block = &(membank->block);
	while (NULL != block) {
		if (block_is_in_buffer(block, membank->type_size, obj)) {
			block_free(block, membank->type_size, obj);
			goto EXIT;
		}
		block = block->next;
	}

	fputs("nb_block: unallocated free.\n", stderr);
	exit(EXIT_FAILURE);
EXIT:
	return;
}

static bool block_is_in_buffer(const block_t *block,
				  uint16_t type_size, char *mem)
{
	return mem >= block->buffer && 
		mem < block->buffer + block->N_max * type_size;

}

static void block_free(block_t *block, uint16_t type_size,
		       char *mem)
{
	int i = (mem - block->buffer) / type_size;
	if (IS_MASK_ENABLED(block->mask, i)) {
		UNSET_MASK(block->mask, i);
		block->N_free += 1;
		set_backward_aligner(block, i);
	} else {
		fputs("nb_block: double free.\n", stderr);
		exit(EXIT_FAILURE);
	}
}

static void set_backward_aligner(block_t *block, int i)
{
	if (i == block->N_align - 1) {
		do {
			block->N_align -= 1;
			if (0 == block->N_align)
				goto EXIT;
		} while (IS_MASK_DISABLED(block->mask,
					  block->N_align - 1));
	}
EXIT:
	return;
}

void nb_membank_finish(nb_membank_t *membank)
{
	block_t *block = membank->block.next;
	while (NULL != block) {
		block_t *next = block->next;
		nb_free(block);
		block = next;
	}
}
