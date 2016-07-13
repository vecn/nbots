#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

#include "nb/memory_bot/nb_malloc.h"
#include "nb/memory_bot/membank.h"

#define MAX(a,b) (((a)>(b))?(a):(b))

#define STATIC_SIZE 8000
#define MASK_SIZE 1000
#define DIV_FREED_TO_ALLOC 3
#define MIN_MAX_DYNAMIC 100

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
static void set_mask(char *mask, int i);
static bool mask_is_set(char *mask, int i);
static void unset_mask(char *mask, int i);
static bool mask_is_unset(char *mask, int i);
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
	while (NULL == mem && NULL != block) {
		mem = block_calloc(block, membank->type_size);
		if (NULL == mem && NULL == block->next)
			block->next = block_create(membank->type_size,
						   membank->N_max_dynamic);
		block = block->next;
	}
	return mem;
}

static void *block_calloc(block_t *block, uint16_t type_size)
{
	int i = block->N_max;
	if (block->N_align < block->N_max) {
		i = block->N_align;
		block->N_align += 1;
	} else if (block->N_free > block->N_max / DIV_FREED_TO_ALLOC) {
		i = 0;
		while (mask_is_set(block->mask, i))
			i += 1;
	}
	void *mem = NULL;
	if (i < block->N_max) {
		mem = block->buffer + i * type_size;
		memset(mem, 0, type_size);
		set_mask(block->mask, i);
		block->N_free -= 1;
	}
	return mem;
}

static block_t *block_create(uint16_t type_size, uint16_t N_max)
{
	uint32_t block_size = sizeof(block_t);
	uint32_t buffer_size = type_size * N_max;
	uint32_t mask_size = N_max / 8 + (N_max % 8 > 0)?1:0;
	uint32_t size = block_size + buffer_size + mask_size;

	char *memory = nb_calloc(size);

	block_t *block = (void*) memory;
	block->N_max = N_max;
	block->N_free = N_max;
	block->buffer = memory + block_size;
	block->mask = memory + block_size + buffer_size;
	return block;
}

static void set_mask(char *mask, int i)
{
	int char_id = i / 8;
	int bit_id = i % 8;
	mask[char_id] = mask[char_id] | (1 << bit_id);
}

static bool mask_is_set(char *mask, int i)
{
	int char_id = i / 8;
	int bit_id = i % 8;
	return (mask[char_id] & (1 << bit_id)) != 0;
}

void nb_membank_free(nb_membank_t *membank, void *obj)
{
	block_t *block = &(membank->block);
	while (NULL != block) {
		if (block_is_in_buffer(block, membank->type_size, obj)) {
			block_free(block, membank->type_size, obj);
			break;
		}
		block = block->next;
	}

	if (NULL == block) {
		fputs("nb_block: unallocated free.\n", stderr);
		exit(EXIT_FAILURE);
	}
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
	if (mask_is_set(block->mask, i)) {
		unset_mask(block->mask, i);
		block->N_free += 1;
		if (i == block->N_align - 1) {
			do {
				block->N_align -= 1;
				if (0 == block->N_align)
					break;
			} while (mask_is_unset(block->mask,
					       block->N_align - 1));
		}
	} else {
		fputs("nb_block: double free.\n", stderr);
		exit(EXIT_FAILURE);
	}
}

static void unset_mask(char *mask, int i)
{
	int char_id = i / 8;
	int bit_id = i % 8;
	mask[char_id] = mask[char_id] & ~(1 << bit_id);
}

static bool mask_is_unset(char *mask, int i)
{
	return !mask_is_set(mask, i);
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
