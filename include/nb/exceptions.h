/**
 * @file exceptions.h
 * @brief Exceptions handlers.
 * @author Victor Eduardo Cardoso Nungaray
 * @n victorc@@cimat.mx
 * @n @@victore_cardoso
 * @date December 8, 2015
 */

#ifndef __NB_EXCEPTIONS_H__
#define __NB_EXCEPTIONS_H__

#include <setjmp.h>

#define nb_exception_try(exception) \
  (0 == setjmp(nb_exception_try_USE_MACRO_INSTEAD_OF_THIS(exception)))

typedef struct nb_exception_s nb_exception_t;

nb_exception_t* nb_exception_create();
void nb_exception_destroy(nb_exception_t *exception);
void nb_exception_throw(nb_exception_t *exception, int id, void *info);
int nb_exception_get_id(nb_exception_t *exception);
char* nb_exception_get_info(nb_exception_t *exception);
void nb_exception_set_alloc(nb_exception_t *exception,
			     void* ptr, void (*free)(void*));
void nb_exception_clear_alloc(nb_exception_t *exception);


void* nb_exception_try_USE_MACRO_INSTEAD_OF_THIS(nb_exception_t *exception);

#endif
