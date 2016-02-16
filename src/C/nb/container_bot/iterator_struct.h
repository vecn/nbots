#ifndef __NB_CONTAINER_BOT_ITERATOR_STRUCT_H__
#define __NB_CONTAINER_BOT_ITERATOR_STRUCT_H__

typedef struct {
	void (*init)(void*);
	void (*copy)(void*, const void*);
	void (*finish)(void*);
	void* (*create)(void);
	void* (*clone)(const void *);
	void (*destroy)(void*);
	void (*clear)(void*);
	void (*restart)(void*);
	const void* (*get_next)(void*);
	bool (*has_more)(const void *);
	void (*set_dst)(void*, const void*);
} iter_binding;

struct nb_iterator_s {
	void *internal;
	iter_binding b;
};

#endif
