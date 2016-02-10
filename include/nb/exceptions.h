/**
 * @file exceptions.h
 * @brief Exceptions handlers.
 * @author Victor Eduardo Cardoso Nungaray
 * @n victorc@@cimat.mx
 * @n @@victore_cardoso
 * @date December 8, 2015
 */

#ifndef __cplusplus /* Do not use with C++ */

#ifndef __VCN_EXCEPTIONS_H__
#define __VCN_EXCEPTIONS_H__

#include <setjmp.h>

#define vcn_exception_try(exception) \
  (0 == setjmp(vcn_exception_try_USE_MACRO_INSTEAD_OF_THIS(exception)))

typedef struct vcn_exception_s vcn_exception_t;

vcn_exception_t* vcn_exception_create();
void vcn_exception_destroy(vcn_exception_t *exception);
void vcn_exception_throw(vcn_exception_t *exception, int id, void *info);
int vcn_exception_get_id(vcn_exception_t *exception);
char* vcn_exception_get_info(vcn_exception_t *exception);
void vcn_exception_set_alloc(vcn_exception_t *exception,
			     void* ptr, void (*free)(void*));
void vcn_exception_clear_alloc(vcn_exception_t *exception);


void* vcn_exception_try_USE_MACRO_INSTEAD_OF_THIS(vcn_exception_t *exception);

#endif

#endif
