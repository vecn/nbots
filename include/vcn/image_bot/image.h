#ifndef __VCN_IMAGE_BOT_IMAGE_H__
#define __VCN_IMAGE_BOT_IMAGE_H__

#include <stdint.h>

typedef struct vcn_image_s vcn_image_t;

vcn_image_t* vcn_image_create(void);
void vcn_image_destroy(vcn_image_t *img);
void vcn_image_read(vcn_image_t *img, const char* filename);
uint32_t vcn_image_get_width(const vcn_image_t *const img);
uint32_t vcn_image_get_height(const vcn_image_t *const img);
uint8_t vcn_image_get_N_channels(const vcn_image_t *const img);
void vcn_image_get_pixel(const vcn_image_t *const img,
			 uint8_t pixel[]);

#endif
