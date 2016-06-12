#ifndef __NB_IMAGE_BOT_IMAGE_H__
#define __NB_IMAGE_BOT_IMAGE_H__

#include <stdint.h>

typedef struct vcn_image_s vcn_image_t;

uint32_t vcn_image_get_memsize(void);
void vcn_image_init(vcn_image_t *img);
void vcn_image_finish(vcn_image_t *img);
void vcn_image_clear(vcn_image_t *img);
vcn_image_t* vcn_image_create(void);
void vcn_image_destroy(vcn_image_t *img);
void vcn_image_init_white(vcn_image_t *img,
			  int width, int height,
			  int comp_x_pixel);
void vcn_image_read(vcn_image_t *img, const char* filename);
uint32_t vcn_image_get_width(const vcn_image_t *const img);
uint32_t vcn_image_get_height(const vcn_image_t *const img);
uint8_t vcn_image_get_N_channels(const vcn_image_t *const img);
void vcn_image_get_pixel(const vcn_image_t *const img, uint32_t r,
			 uint32_t c,  uint8_t pixel[]);
void vcn_image_resize(const vcn_image_t *input_img,
		      vcn_image_t *output_img,
		      int out_width, int out_height);
void vcn_image_write_png(const vcn_image_t *img, const char *filename);
void vcn_image_render_ascii(const vcn_image_t *img, const char *filename,
			    int pixels_x_char);

#endif
