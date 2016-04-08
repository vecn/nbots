#ifndef __NB_CONTAINER_BOT_CONTAINER_STRUCT_H__
#define __NB_CONTAINER_BOT_CONTAINER_STRUCT_H__

#include <stdint.h>
#include <stdbool.h>

#include "nb/container_bot/container.h"

typedef struct {
	uint32_t (*key)(const void*); 
	void (*destroy)(void*);
	int8_t (*compare)(const void*, const void*);
	void* (*clone)(const void*);
} value_operators;

typedef struct {
	void (*init)(void*);
	void (*copy)(void*, const void*,
		     void* (*clone)(const void*));
	void (*finish)(void*,
		       void (*destroy)(void*));
	void* (*create)(void);
	void* (*clone)(const void *,
		       void* (*clone)(const void*));
	void (*destroy)(void*,
			void (*destroy)(void*));
	void (*clear)(void*,
		      void (*destroy)(void*));
	void (*merge)(void*, void*,
		      uint32_t (*key)(const void*),
		      int8_t (*compare)(const void*, const void*));
	bool (*insert)(void*, const void *const,
		       uint32_t (*key)(const void*),
		       int8_t (*compare)(const void*, const void*));
	void* (*get_first)(const void* const);
	void* (*delete_first)(void*,
			      uint32_t (*key)(const void*));
	void* (*exist)(const void*, const void*,
		       uint32_t (*key)(const void*),
		       int8_t (*compare)(const void*, const void*));
	void* (*delete)(void*, const void*,
			uint32_t (*key)(const void*),
			int8_t (*compare)(const void*, const void*));
	uint32_t (*get_length)(const void*);
	bool (*is_empty)(const void*);
	bool (*is_not_empty)(const void*);	
} container_binding;

struct nb_container_s {
	nb_container_type type;
	void *cnt;
	value_operators op;
	container_binding b;
};

#endif
