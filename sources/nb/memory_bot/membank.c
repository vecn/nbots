#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <signal.h>

#include "nb/memory_bot/allocate_mem.h"
#include "nb/memory_bot/membank.h"

#define DEFAULT_N_MAX 1024

#define GET_64_ALIGNMENT(N_max)			\
	((N_max)/64 + (((N_max)%64 > 0)?1:0))
#define GET_MASK_SIZE(N_max)			\
	(8 * GET_64_ALIGNMENT((N_max)))
#define SET_MASK(mask, i)			\
	((mask)[(i)/8] |= 1 << ((i)%8))
#define IS_MASK_ENABLED(mask, i)		\
	((mask)[(i)/8] & (1 << ((i)%8)))
#define UNSET_MASK(mask, i)			\
	((mask)[(i)/8] &=  ~(1 << ((i)%8)))
#define IS_MASK_DISABLED(mask, i)		\
	(!IS_MASK_ENABLED(mask, i))
#define IS_MASK_64_FULL(mask, i)		\
	(~((uint64_t*)(mask))[(i)] == 0)
#define IS_MASK_32_FULL(mask, i)		\
	(~((uint32_t*)(mask))[(i)] == 0)
#define IS_MASK_16_FULL(mask, i)		\
	(~((uint16_t*)(mask))[(i)] == 0)
#define IS_MASK_BYTE_FULL(mask, i)		\
	(~((mask)[(i)]) == 0)
#define IS_MASK_64_EMPTY(mask, i)		\
	(((uint64_t*)(mask))[(i)] == 0)
#define IS_MASK_32_EMPTY(mask, i)		\
	(((uint32_t*)(mask))[(i)] == 0)
#define IS_MASK_16_EMPTY(mask, i)		\
	(((uint16_t*)(mask))[(i)] == 0)
#define IS_MASK_BYTE_EMPTY(mask, i)		\
	((mask)[(i)] == 0)

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
	block_t *block;
	uint16_t type_size;
	uint16_t N_max;
};

static void* block_data_calloc(block_t *block, uint16_t type_size);
static int get_first_free_id(const char *mask);
static block_t *block_create(uint16_t type_size, uint16_t N_max);
static bool block_is_in_buffer(const block_t *block,
			       uint16_t type_size, char *mem);
static void block_data_free(block_t *block, uint16_t type_size,
			    char *mem);
static int get_last_allocated_id(const char *mask, uint16_t N_max);
static void block_clear(block_t *block);

uint32_t nb_membank_get_memsize(void)
{
	return sizeof(nb_membank_t);
}

void nb_membank_init(nb_membank_t *membank, uint16_t type_size)
{
	memset(membank, 0, sizeof(*membank));
	membank->type_size = type_size;
	membank->N_max = DEFAULT_N_MAX;
}

void *nb_membank_allocate_mem(nb_membank_t *membank)
{
	block_t *block;
	if (NULL == membank->block) {
		block = block_create(membank->type_size, membank->N_max);
		membank->block = block;
	} else {
		block = membank->block;
		while (0 == block->N_free) {
			if (NULL == block->next)
				block->next = block_create(membank->type_size,
							   membank->N_max);
			block = block->next;
		}
	}
	return block_data_calloc(block, membank->type_size);
}

static void *block_data_calloc(block_t *block, uint16_t type_size)
{
	void *mem = NULL;
	int i;
	if (block->N_align < block->N_max) {
		i = block->N_align;
		block->N_align += 1;
	} else if (block->N_free > 0) {
		i = get_first_free_id(block->mask);
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

static int get_first_free_id(const char *mask)
{
	int i = 0;
	while (IS_MASK_64_FULL(mask, i))
		i++;

	i <<= 1;
	while (IS_MASK_32_FULL(mask, i))
		i++;

	i <<= 1;
	while (IS_MASK_16_FULL(mask, i))
		i++;

	i <<= 1;
	while (IS_MASK_BYTE_FULL(mask, i))
		i++;

	i <<= 3;
	while (IS_MASK_ENABLED(mask, i))
		i++;
	return i;
}

static block_t *block_create(uint16_t type_size, uint16_t N_max)
{
	uint32_t block_size = sizeof(block_t);
	uint32_t buffer_size = type_size * N_max;
	uint32_t mask_size = GET_MASK_SIZE(N_max);
	uint32_t size = block_size + buffer_size + mask_size;

	char *memory = nb_allocate_zero_mem(size);

	block_t *block = (void*) memory;
	block->N_max = N_max;
	block->N_free = N_max;
	block->buffer = memory + block_size;
	block->mask = memory + block_size + buffer_size;
	return block;
}


void nb_membank_free_mem(nb_membank_t *membank, void *obj)
{
	block_t *block = membank->block;
	while (NULL != block) {
		if (block_is_in_buffer(block, membank->type_size, obj)) {
			block_data_free(block, membank->type_size, obj);
			goto EXIT;
		}
		block = block->next;
	}

	fputs("\nNB_MEMORY ERROR: unallocated free.\n\n", stderr);
	raise(SIGSEGV);
EXIT:
	return;
}

static bool block_is_in_buffer(const block_t *block,
			       uint16_t type_size, char *mem)
{
	return mem >= block->buffer && 
		mem < block->buffer + block->N_max * type_size;

}

static void block_data_free(block_t *block, uint16_t type_size,
			    char *mem)
{
	int i = (mem - block->buffer) / type_size;
	if (IS_MASK_ENABLED(block->mask, i)) {
		UNSET_MASK(block->mask, i);
		block->N_free += 1;
		if (i == block->N_align - 1)
			block->N_align =
				1 + get_last_allocated_id(block->mask,
							  block->N_max);
	} else {
		fputs("\nNB_MEMORY ERROR: double free.\n\n", stderr);
		raise(SIGSEGV);
	}
}

static int get_last_allocated_id(const char *mask, uint16_t N_max)
{
	int i = GET_64_ALIGNMENT(N_max) - 1;
	while (IS_MASK_64_EMPTY(mask, i)) {
		i--;
		if (i < 0)
			goto EXIT;
	}

	i = (i << 1) + 1;
	while (IS_MASK_32_EMPTY(mask, i))
		i--;

	i = (i << 1) + 1;
	while (IS_MASK_16_EMPTY(mask, i))
		i--;

	i = (i << 1) + 1;
	while (IS_MASK_BYTE_EMPTY(mask, i))
		i--;

	i = (i << 3) + 7;
	while (IS_MASK_DISABLED(mask, i))
		i--;
EXIT:
	return i;
}

void nb_membank_finish(nb_membank_t *membank)
{
	block_t *block = membank->block;
	while (NULL != block) {
		block_t *next = block->next;
		nb_free_mem(block);
		block = next;
	}
}

void nb_membank_clear(nb_membank_t *membank)
{
	block_t *block = membank->block;
	while (NULL != block) {
		block_clear(block);
		block = block->next;
	}
}

static void block_clear(block_t *block)
{
	block->N_align = 0;
	block->N_free = block->N_max;
	uint32_t mask_size = GET_MASK_SIZE(block->N_max);
	memset(block->mask, 0, mask_size);
}

void nb_membank_merge(nb_membank_t *membank, nb_membank_t *extension)
{
	if (NULL != membank->block) {
		block_t *block = membank->block;
		while (NULL != block->next)
			block = block->next;
		block->next = extension->block;
	} else {
		membank->block = extension->block;
	}
	extension->block = NULL;
}

void nb_membank_set_N_x_block(nb_membank_t *membank, uint16_t N_x_block)
{
	membank->N_max = 64 * GET_64_ALIGNMENT(N_x_block);
}
