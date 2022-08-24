#include "nb/memory_bot.h"
#include "nb/container_bot.h"

#include "dewall_dependencies.h"

static uint32_t hash_table_size(void)
{
	return nb_container_get_memsize(NB_HASH);
}

static void hash_table_init(afl_t* self,
			    uint32_t (*keygen)(const void *const))
{
	nb_container_init((nb_container_t*)self, NB_HASH);
	nb_container_set_key_generator((nb_container_t*)self, keygen);
}

static void hash_table_finish(afl_t* self)
{
	nb_container_finish((nb_container_t*)self);
}

static bool hash_table_is_empty(const afl_t *const self)
{
	return nb_container_is_empty((const nb_container_t *const)self);
}

static void hash_table_insert(afl_t* self, const void *const elem)
{
	nb_container_insert((nb_container_t*)self, elem);
}

static void* hash_table_delete_any(afl_t* self)
{
	return nb_container_delete_first((nb_container_t*)self);
}

static void* hash_table_delete(afl_t* self, const void *const elem)
{
	return nb_container_delete((nb_container_t*)self, elem);
}

static uint32_t hash_table_iterator_size(void)
{
	return nb_iterator_get_memsize();
}

static void hash_table_iterator_init(afl_iterator_t* iter)
{
	nb_iterator_init((nb_iterator_t*)iter);
}

static void hash_table_iterator_finish(afl_iterator_t* iter)
{
	nb_iterator_finish((nb_iterator_t*)iter);
}

static void hash_table_iterator_set_afl(afl_iterator_t* iter,
					const afl_t *const afl)
{
	nb_iterator_set_container((nb_iterator_t*)iter,
				  (const nb_container_t *const)afl);
}

static const void* hash_table_iterator_get_next(afl_iterator_t* iter)
{
	return nb_iterator_get_next((nb_iterator_t*)iter);
}

static bool hash_table_iterator_has_more(const afl_iterator_t *const iter)
{
	return nb_iterator_has_more((nb_iterator_t*)iter);
}

static uint32_t queue_size(void)
{
	return nb_container_get_memsize(NB_QUEUE);
}

static void queue_init(queue_t* self)
{
	nb_container_init((nb_container_t*)self, NB_QUEUE);
}

static void queue_finish(queue_t* self)
{
	nb_container_finish((nb_container_t*)self);
}

static void queue_add(queue_t *self, const void *const elem)
{
	nb_container_insert((nb_container_t*)self, elem);
}

static void* queue_poll(queue_t *self)
{
	return nb_container_delete_first((nb_container_t*)self);
}

static bool queue_is_empty(const queue_t *const self)
{
	return nb_container_is_empty((nb_container_t*)self);
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
		hash_table_is_empty,
		hash_table_insert,
		hash_table_delete_any,
		hash_table_delete
	},
	{
		hash_table_iterator_size,
		hash_table_iterator_init,
		hash_table_iterator_finish,
		hash_table_iterator_set_afl,
		hash_table_iterator_get_next,
		hash_table_iterator_has_more,
	},
	{
		queue_size,
		queue_init,
		queue_finish,
		queue_add,
		queue_poll,
		queue_is_empty
	}
};

interface_t* module(void)
{
	return &implementation;
}
