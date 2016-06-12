#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "nb/image_bot/image.h"

#define STB_IMAGE_IMPLEMENTATION
#include "tiny_libs/stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "tiny_libs/stb_image_write.h"

#define STB_IMAGE_RESIZE_IMPLEMENTATION
#include "tiny_libs/stb_image_resize.h"

struct vcn_image_s {
	uint32_t width;
	uint32_t height;
	uint8_t comp_x_pixel;
	uint8_t* pixels;
};

uint32_t vcn_image_get_memsize(void)
{
	return sizeof(vcn_image_t);
}

void vcn_image_init(vcn_image_t *img)
{
	uint32_t memsize = vcn_image_get_memsize();
	memset(img, 0, memsize);
}

void vcn_image_finish(vcn_image_t *img)
{
	vcn_image_clear(img);
}

void vcn_image_clear(vcn_image_t *img)
{
	if (NULL != img->pixels)
		stbi_image_free(img->pixels);
	uint32_t memsize = vcn_image_get_memsize();
	memset(img, 0, memsize);
}

inline vcn_image_t* vcn_image_create(void)
{
	uint32_t memsize = vcn_image_get_memsize();
	return calloc(1, memsize);
}

void vcn_image_destroy(vcn_image_t *img)
{
	vcn_image_clear(img);
	free(img);
}

void vcn_image_init_white(vcn_image_t *img,
			  int width, int height,
			  int comp_x_pixel)
{
	img->width = width;
	img->height = height;
	img->comp_x_pixel;
	uint32_t memsize = width * height * comp_x_pixel;
	img->pixels = calloc(memsize, 1);
}

void vcn_image_read(vcn_image_t *img, const char* filename)
{
	int w, h, n;
	uint8_t *const restrict pixels = stbi_load(filename, &w, &h, &n, 0);

	if (pixels != NULL) {
		img->pixels = pixels;
		img->width = w;
		img->height = h;
		img->comp_x_pixel = n;
	}
}

inline uint32_t vcn_image_get_width(const vcn_image_t *const img)
{
	return img->width;
}

inline uint32_t vcn_image_get_height(const vcn_image_t *const img)
{
	return img->height;
}

inline uint8_t vcn_image_get_N_channels(const vcn_image_t *const img)
{
	return img->comp_x_pixel;
}

void vcn_image_get_pixel(const vcn_image_t *const img, uint32_t r,
			 uint32_t c, uint8_t pixel[])
{
	if (NULL != img->pixels) {
		uint32_t w = img->width * img->comp_x_pixel;
		uint32_t col = c * img->comp_x_pixel;
		memcpy(pixel, &(img->pixels[r * w + col]), img->comp_x_pixel);
	}
}

void vcn_image_resize(const vcn_image_t *input_img,
		      vcn_image_t *output_img,
		      int out_width, int out_height)
{
	output_img->width = out_width;
	output_img->height = out_height;
	output_img->comp_x_pixel = input_img->comp_x_pixel;
	uint32_t memsize = input_img->comp_x_pixel * out_width * out_height;
	output_img->pixels = malloc(memsize);

	int stride_in_bytes = input_img->width * input_img->comp_x_pixel;
	int output_stride_in_bytes = out_width * input_img->comp_x_pixel;
	int status = stbir_resize_uint8(input_img->pixels,
					input_img->width,
					input_img->height,
					stride_in_bytes,
					output_img->pixels,
					out_width, out_height,
					output_stride_in_bytes,
					input_img->comp_x_pixel);
	assert(0 != status);
}

void vcn_image_write_png(const vcn_image_t *img, const char *filename)
{  
	int stride_in_bytes = img->width * img->comp_x_pixel;
	int status = stbi_write_png(filename,
				    img->width,
				    img->height,
				    img->comp_x_pixel,
				    img->pixels,
				    stride_in_bytes);
	assert(0 != status);
}

void vcn_image_render_ascii(const vcn_image_t *img, const char *filename,
			    int pixels_x_char)
{
	;
	/* LUMA conversion: TEMPORAL */
  
	/* Long gray scale */
	//"$@B%8&WM#*";
	//"oahkbdpqwm";
	//"ZO0QLCJUYX";
	//"zcvunxrjft";
	//"/\|()1{}[]";
	//"?-_+~<>i!l";
	//"I;:,\"^`'. ";

	/* Short gray scale */
	//"@%#*+=:-. ";

	/* Faltan: */
	//"egsyADEFGH";
	//"KNPRSTV234";
	//"5679=";
}
