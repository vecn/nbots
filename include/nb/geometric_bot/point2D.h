#ifndef __NB_GEOMETRIC_BOT_POINT2D_H__
#define __NB_GEOMETRIC_BOT_POINT2D_H__

#include <stdint.h>
#include <stdbool.h>
#include "nb/container_bot/container.h"

typedef struct vcn_point2D_s {
	double x[2];
	void* attr;
} vcn_point2D_t;

vcn_point2D_t*  vcn_point2D_create(void);
void vcn_point2D_destroy(void *point_ptr);
int8_t vcn_point2D_compare(const void *const p1_ptr,
			   const void *const p2_ptr);

#endif
