#ifndef __VCN_GEOMETRIC_BOT_POINT2D_H__
#define __VCN_GEOMETRIC_BOT_POINT2D_H__

#include <stdint.h>
#include <stdbool.h>
#include "vcn/container_bot/container.h"

#ifdef __cplusplus
extern "C" {
#endif
	typedef struct vcn_point2D_s {
		double x[2];
		void* attr;
	} vcn_point2D_t;

	vcn_point2D_t*  vcn_point2D_create(void);
	void vcn_point2D_destroy(vcn_point2D_t* point);
	bool vcn_point2D_are_equal(const void *const p1_ptr,
				   const void *const p2_ptr);

#ifdef __cplusplus
}
#endif

#endif
