#include "nb/memory_bot.h"
#include "nb/container_bot.h"

#include "dewall_dependencies.h"

static uint32_t hash_table_size(void)
{
	return nb_container_get_memsize(NB_HASH);
}

static void hash_table_init(hash_table_t* self)
{
	nb_container_init((nb_container_t*)self, NB_HASH);
}

static void hash_table_finish(hash_table_t* self)
{
	nb_container_finish((nb_container_t*)self);
}

static void hash_table_set_keygen(hash_table_t* self,
				  uint32_t (*keygen)(const void *const))
{
	nb_container_set_key_generator((nb_container_t*)self, keygen);
}

static bool hash_table_is_empty(const hash_table_t *const self)
{
	return nb_container_is_empty((const nb_container_t *const)self);
}

static void hash_table_insert(hash_table_t* self, const void *const elem)
{
	nb_container_insert((nb_container_t*)self, elem);
}

static void* hash_table_delete_any(hash_table_t* self)
{
	return nb_container_delete_first((nb_container_t*)self);
}

static void* hash_table_delete(hash_table_t* self, const void *const elem)
{
	return nb_container_delete((nb_container_t*)self, elem);
}

static interface_t implementation = {
	{
		nb_allocate_mem,
		nb_allocate_zero_mem,
		nb_free_mem
	},
	{
		hash_table_size,
		hash_table_init,
		hash_table_finish,
		hash_table_set_keygen,
		hash_table_is_empty,
		hash_table_insert,
		hash_table_delete_any,
		hash_table_delete
	}
};

interface_t* module(void)
{
	return &implementation;
}
