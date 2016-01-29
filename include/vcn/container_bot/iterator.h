/**
 * @file iterator.h
 * @brief  Generic iterator.
 * @author Victor Eduardo Cardoso Nungaray
 * @n victorc@@cimat.mx
 * @n @@victore_cardoso
 */

#ifndef __VCN_ITERATOR_H__
#define __VCN_ITERATOR_H__


#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include "vcn/container_bot/container.h"
  
	typedef struct vcn_iterator_s vcn_iterator_t;

	void vcn_iterator_init(vcn_iterator_t *iter);
	vcn_iterator_t* vcn_iterator_create(void);
	void vcn_iterator_set_container(vcn_iterator_t *iter,
					const vcn_container_t *const container);
	vcn_iterator_t* vcn_iterator_clone(const vcn_iterator_t *const iter);
	void vcn_iterator_destroy(vcn_iterator_t *iter);
	void vcn_iterator_restart(vcn_iterator_t *iter);
	
	const void* vcn_iterator_get_next(vcn_iterator_t *iter);
	bool vcn_iterator_has_more(const vcn_iterator_t *const iter);

#ifdef __cplusplus
}
#endif

#endif
