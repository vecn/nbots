#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/geometric_bot/mesh/tessellator2D.h"

#include "tessellator2D_structs.h"

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

static void hash_table_iterator_init(afl_iterator_t* iter,
				     const afl_t *const afl)
{
	nb_iterator_init((nb_iterator_t*)iter);
	nb_iterator_set_container((nb_iterator_t*)iter,
				  (const nb_container_t *const)afl);
}

static void hash_table_iterator_finish(afl_iterator_t* iter)
{
	nb_iterator_finish((nb_iterator_t*)iter);
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

static uint32_t mesh_edge_iterator_size(void)
{
	return nb_iterator_get_memsize();
}

static void mesh_edge_iterator_init(edge_iterator_t* iter,
				    const mesh_t *const mesh)
{
	nb_tessellator2D_t* m2D = (nb_tessellator2D_t*) mesh;
	nb_iterator_init((nb_iterator_t*)iter);
	nb_iterator_set_container((nb_iterator_t*)iter, m2D->ht_edge);
}

static void mesh_edge_iterator_finish(edge_iterator_t* iter)
{
	nb_iterator_finish((nb_iterator_t*)iter);
}

static const void* mesh_edge_iterator_get_next(
	const edge_iterator_t *const iter)
{
	return nb_iterator_get_next((nb_iterator_t*)iter);
}

static bool mesh_edge_iterator_has_more(const edge_iterator_t *const iter)
{
	return nb_iterator_has_more((nb_iterator_t*)iter);
}

static trg_t* mesh_new_triangle(mesh_t* self, /* Counter-clockwise order */
				const vtx_t *const v1,
				const vtx_t *const v2,
				const vtx_t *const v3)
{
	nb_tessellator2D_t *mesh = (nb_tessellator2D_t*) self;
	msh_trg_t *trg = mtrg_allocate_zero_mem(mesh);
	trg->v1 = (msh_vtx_t*) v1;
	trg->v2 = (msh_vtx_t*) v2;
	trg->v3 = (msh_vtx_t*) v3;
	return (trg_t*) trg;
}

static void mesh_connect_triangle(mesh_t *self, const trg_t *const trg)
{
	mesh_add_triangle((nb_tessellator2D_t*)self, (msh_trg_t*)trg);
}

static void on_triangle_connection(const mesh_t *const self,
				   const trg_t *const trg)
{
	/* Useful to create animations */
	nb_tessellator2D_t *mesh = (nb_tessellator2D_t*) self;
	mesh->do_after_insert_trg(mesh);
}

// TODO: Separate static structure
static interface_t implementation = {
	{
		/* interface_heap_t mem; */
		nb_allocate_mem,
		nb_allocate_zero_mem,
		nb_free_mem
	},
	{
		/* interface_afl_t afl */
		hash_table_size,
		hash_table_init,
		hash_table_finish,
		hash_table_is_empty,
		hash_table_insert,
		hash_table_delete_any,
		hash_table_delete
	},
	{
		/* interface_afl_iterator_t afl_iter */
		hash_table_iterator_size,
		hash_table_iterator_init,
		hash_table_iterator_finish,
		hash_table_iterator_get_next,
		hash_table_iterator_has_more,
	},
	{
		/* interface_queue_t queue */
		queue_size,
		queue_init,
		queue_finish,
		queue_add,
		queue_poll,
		queue_is_empty
	},
	{
		/* interface_mesh_t mesh */
		mesh_new_triangle,
		mesh_connect_triangle,
		on_triangle_connection,
		{
			/* interface_mesh_edge_iterator_t edge_iter */
			mesh_edge_iterator_size,
			mesh_edge_iterator_init,
			mesh_edge_iterator_finish,
			mesh_edge_iterator_get_next,
			mesh_edge_iterator_has_more
		}
	}
};

interface_t* module(void)
{
	return &implementation;
}
