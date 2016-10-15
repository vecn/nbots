#ifndef __NB_IMAGE_BOT_IMAGE_H__
#define __NB_IMAGE_BOT_IMAGE_H__

#include <stdint.h>

typedef struct {
	uint32_t width;
	uint32_t height;
	uint8_t comp_x_pixel;
	/* 1: grey
	 * 2: grey, alpha
	 * 3: red, green, blue
	 * 4: red, green, blue, alpha
	 */
	uint8_t* pixels;
} nb_image_t;

uint32_t nb_image_get_memsize(void);
void nb_image_init(nb_image_t *img);
void nb_image_finish(nb_image_t *img);
void nb_image_clear(nb_image_t *img);
void nb_image_init_white(nb_image_t *img,
			 int width, int height,
			 int comp_x_pixel);
void nb_image_read(nb_image_t *img, const char* filename);
uint32_t nb_image_get_width(const nb_image_t *const img);
uint32_t nb_image_get_height(const nb_image_t *const img);
uint8_t nb_image_get_N_channels(const nb_image_t *const img);
void nb_image_get_pixel(const nb_image_t *const img, uint32_t r,
			uint32_t c,  uint8_t pixel[]);
void nb_image_set_pixel(nb_image_t *img, uint32_t r,
			uint32_t c,  uint8_t pixel[]);
void nb_image_blend_pixel_ga(nb_image_t *img, uint32_t r,
			     uint32_t c,  uint8_t pixel[2]);
void nb_image_blend_pixel_rgba(nb_image_t *img, uint32_t r,
			       uint32_t c,  uint8_t pixel[4]);
uint8_t nb_image_get_pixel_luma(const nb_image_t *img,
				uint32_t r, uint32_t c);
void nb_image_resize(const nb_image_t *input_img,
		     nb_image_t *output_img,
		     int out_width, int out_height);
void nb_image_write(const nb_image_t *img, const char *filename);
void nb_image_write_png(const nb_image_t *img, const char *filename);
void nb_image_write_bmp(const nb_image_t *img, const char *filename);
void nb_image_write_tga(const nb_image_t *img, const char *filename);
void nb_image_write_ascii(const nb_image_t *img, const char *filename,
			  uint16_t N_chars_width);

#endif
