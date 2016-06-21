#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <alloca.h>

#include "nb/image_bot/image.h"

#define STB_IMAGE_STATIC
#define STB_IMAGE_IMPLEMENTATION
#include "tiny_libs/stb_image.h"

#define STB_IMAGE_WRITE_STATIC
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "tiny_libs/stb_image_write.h"

#define STB_IMAGE_RESIZE_STATIC
#define STB_IMAGE_RESIZE_IMPLEMENTATION
#include "tiny_libs/stb_image_resize.h"

#define _LUMA_R_WEIGHT 0.299
#define _LUMA_G_WEIGHT 0.587
#define _LUMA_B_WEIGHT 0.114
#define _ASCII_CHAR_ASPECT_RATIO 0.5

#define GET_FILENAME_EXT(filename, ext)				\
	do {							\
		const char *dot = strrchr((filename), '.');	\
		if (!dot || dot == (filename))			\
			(ext) = "";				\
		else						\
			(ext) = dot + 1;			\
	} while(0)

#define TO_LOWERCASE(str)					\
	do {							\
	        int i = 0;					\
		while ((str)[i] != '\0' && (str)[i] != '\n') {	\
			(str)[i] = tolower((str)[i]);		\
			i++;					\
		}						\
	} while(0)

enum {
	PNG, BMP, TGA, HDR, TXT, UNSUPPORTED
};

static int get_format(const char *filename);
static void image_write_default(const vcn_image_t *img, const char *filename);

static const char *get_ascii_scale(void);
static uint8_t get_block_luma(const vcn_image_t *img, uint32_t i, uint32_t j,
			      uint16_t block_width, uint16_t block_height);
static char get_ascii_from_luma(uint8_t luma, const char *scale,
				uint8_t *next_char);
static uint8_t get_luma_1comp(const vcn_image_t *img, uint32_t r, uint32_t c);
static uint8_t get_luma_2comp(const vcn_image_t *img, uint32_t r, uint32_t c);
static uint8_t get_luma_3comp(const vcn_image_t *img, uint32_t r, uint32_t c);
static uint8_t get_luma_4comp(const vcn_image_t *img, uint32_t r, uint32_t c);

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
	vcn_image_finish(img);
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
	uint32_t w = img->width * img->comp_x_pixel;
	uint32_t col = c * img->comp_x_pixel;
	memcpy(pixel, &(img->pixels[r * w + col]), img->comp_x_pixel);
}

void vcn_image_set_pixel(vcn_image_t *img, uint32_t r,
			 uint32_t c,  uint8_t pixel[])
{
	uint32_t w = img->width * img->comp_x_pixel;
	uint32_t col = c * img->comp_x_pixel;
	memcpy(&(img->pixels[r * w + col]), pixel, img->comp_x_pixel);
}

void vcn_image_blend_pixel_ga(vcn_image_t *img, uint32_t r,
				uint32_t c,  uint8_t pixel[2])
{
	uint32_t w = img->width * img->comp_x_pixel;
	uint32_t col = c * img->comp_x_pixel;
	float dest_alpha = img->pixels[r * w + col + 1] / 255.0f;
	float src_alpha = pixel[1] / 255.0f;
	/* Alpha channel */
	float out_alpha = src_alpha + (1.0 - src_alpha) * dest_alpha;
	img->pixels[r * w + col + 1] = (uint8_t)(255.0 * out_alpha + 0.5);
	/* Grayscale */
	if (img->pixels[r * w + col + 1] > 0) {
		float a = src_alpha / out_alpha;
		float b = dest_alpha * (1.0f - src_alpha) / out_alpha;
		img->pixels[r * w + col] = (uint8_t)
			(a * pixel[0] + b * img->pixels[r * w + col] + 0.5);
	}
}

void vcn_image_blend_pixel_rgba(vcn_image_t *img, uint32_t r,
				uint32_t c,  uint8_t pixel[4])
{
	uint32_t w = img->width * img->comp_x_pixel;
	uint32_t col = c * img->comp_x_pixel;
	float dest_alpha = img->pixels[r * w + col + 3] / 255.0f;
	float src_alpha = pixel[3] / 255.0f;
	/* Alpha channel */
	float out_alpha = src_alpha + (1.0 - src_alpha) * dest_alpha;
	img->pixels[r * w + col + 3] = (uint8_t)(255.0 * out_alpha + 0.5);
	/* RGB */
	if (img->pixels[r * w + col + 3] > 0) {
		float a = src_alpha / out_alpha;
		float b = dest_alpha * (1.0f - src_alpha) / out_alpha;
		for (uint8_t i = 0; i < 3; i++)
			img->pixels[r * w + col + i] = (uint8_t)
				(a * pixel[i] + 
				 b * img->pixels[r * w + col + i] + 0.5);
	}
}

uint8_t vcn_image_get_pixel_luma(const vcn_image_t *img,
				 uint32_t r, uint32_t c)
{
	uint8_t luma;
	switch (img->comp_x_pixel) {
	case 1:
		luma = get_luma_1comp(img, r, c);
		break;
	case 2:
		luma = get_luma_2comp(img, r, c);
		break;
	case 3:
		luma = get_luma_3comp(img, r, c);
		break;
	case 4:
		luma = get_luma_4comp(img, r, c);
		break;
	}
	return  luma;
}

static uint8_t get_luma_1comp(const vcn_image_t *img, uint32_t r, uint32_t c)
{
	uint32_t w = img->width;
	return img->pixels[r * w + c];
}

static uint8_t get_luma_2comp(const vcn_image_t *img, uint32_t r, uint32_t c)
{
	uint32_t w = 2 * img->width;
	uint8_t g = img->pixels[r * w + (c * 2)];
	uint8_t a = img->pixels[r * w + (c*2+1)];
	double weight = a / 255.0;
	return (1 - weight) * 255 + weight * g;
}

static uint8_t get_luma_3comp(const vcn_image_t *img, uint32_t r, uint32_t c)
{
	uint32_t w = 3 * img->width;
	uint8_t red = img->pixels[r * w + (c * 3)];
	uint8_t g = img->pixels[r * w + (c*3+1)];
	uint8_t b = img->pixels[r * w + (c*3+2)];
	return (uint8_t)(_LUMA_R_WEIGHT * red + 0.5) +
		(uint8_t)(_LUMA_G_WEIGHT * g + 0.5) +
		(uint8_t)(_LUMA_B_WEIGHT * b + 0.5);
}

static uint8_t get_luma_4comp(const vcn_image_t *img, uint32_t r, uint32_t c)
{
	uint32_t w = 4 * img->width;
	uint8_t red = img->pixels[r * w + (c * 4)];
	uint8_t g = img->pixels[r * w + (c*4+1)];
	uint8_t b = img->pixels[r * w + (c*4+2)];
	uint8_t a = img->pixels[r * w + (c*4+3)];
	double weight = a / 255.0;
	return (1 - weight) * 255 +
		weight * ((uint8_t)(_LUMA_R_WEIGHT * red + 0.5) +
			  (uint8_t)(_LUMA_G_WEIGHT * g + 0.5) +
			  (uint8_t)(_LUMA_B_WEIGHT * b + 0.5));
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

void vcn_image_write(const vcn_image_t *img, const char *filename)
{
	int format = get_format(filename);
	switch (format) {
	case PNG:
		vcn_image_write_png(img, filename);
		break;
	case BMP:
		vcn_image_write_bmp(img, filename);
		break;
	case TGA:
		vcn_image_write_tga(img, filename);
		break;
	case HDR:
		vcn_image_write_hdr(img, filename);
		break;
	case TXT:
		vcn_image_write_ascii(img, filename,
				      img->width / 10);
		break;
	default:
		image_write_default(img, filename);
	}
}

static int get_format(const char *filename)
{
	const char *raw_ext;
	GET_FILENAME_EXT(filename, raw_ext);
	char ext[10];
	strcpy(ext, raw_ext);
	TO_LOWERCASE(ext);
	int format;
	if (0 == strcmp(ext, "png"))
		format = PNG;
	else if (0 == strcmp(ext, "bmp"))
		format = BMP;
	else if (0 == strcmp(ext, "tga"))
		format = TGA;
	else if (0 == strcmp(ext, "hdr"))
		format = HDR;
	else if (0 == strcmp(ext, "txt"))
		format = TXT;
	else
		format = UNSUPPORTED;
	return format;
}

static void image_write_default(const vcn_image_t *img, const char *filename)
{
	char *name = alloca(strlen(filename) + 4);
	sprintf(name, "%s.png", filename);
	vcn_image_write_png(img, name);
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

void vcn_image_write_bmp(const vcn_image_t *img, const char *filename)
{  
	int status = stbi_write_bmp(filename,
				    img->width,
				    img->height,
				    img->comp_x_pixel,
				    img->pixels);
	assert(0 != status);
}

void vcn_image_write_hdr(const vcn_image_t *img, const char *filename)
{  
	int status = stbi_write_hdr(filename,
				    img->width,
				    img->height,
				    img->comp_x_pixel,
				    img->pixels);
	assert(0 != status);
}

void vcn_image_write_tga(const vcn_image_t *img, const char *filename)
{  
	int status = stbi_write_tga(filename,
				    img->width,
				    img->height,
				    img->comp_x_pixel,
				    img->pixels);
	assert(0 != status);
}

void vcn_image_write_ascii(const vcn_image_t *img, const char *filename,
			   uint16_t N_chars_width)
{
	FILE *fp = fopen(filename, "w");
	if (NULL == fp)
		goto EXIT;

	const char *scale = get_ascii_scale();
	uint16_t block_width = img->width / N_chars_width;
	uint16_t block_height = block_width / _ASCII_CHAR_ASPECT_RATIO;
	uint32_t w = img->width / block_width;
	uint32_t h = img->height / block_height;
	uint8_t next_char = 0;
	for (uint32_t i = 0; i < h; i++) {
		for (uint32_t j = 0; j < w; j++) {
			uint8_t luma = get_block_luma(img, i, j,
						      block_width,
						      block_height);
			char c = get_ascii_from_luma(luma, scale,
						     &next_char);
			fprintf(fp, "%c", c);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
EXIT:
	return;
}

static const char *get_ascii_scale(void)
{
	return  "W%o*+~:-. "  \
		"Wgeu+~:-. "  \
		"W$rz+~:-. "  \
		"@#xj+~:-. "  \
		"@#Qv+~:-. "  \
		"B$pc+~:-. "  \
		"Bkqw+~:-. "  \
		"Mhos+~:-. "  \
		"8&bai~:-. ";

}

static uint8_t get_block_luma(const vcn_image_t *img, uint32_t i, uint32_t j,
			      uint16_t block_width, uint16_t block_height)
{
	uint16_t pix_jump = 1;
	uint16_t Nbh = block_height;
	uint16_t Nbw = block_width;
	if (block_height > 10) {
		pix_jump = block_height / 10;
		Nbh = 10;
		Nbw = block_width / pix_jump;
	}
		
	uint32_t sum_luma = 0;	
	for (uint16_t m = 0; m < Nbh; m++) {
		uint32_t r = i * block_height + pix_jump * m;
		for (uint16_t n = 0; n < Nbw; n++) {
			uint32_t c = j * block_width + pix_jump * n;
			sum_luma += vcn_image_get_pixel_luma(img, r, c);			
		}
	}
	return sum_luma / (Nbw * Nbh);
}

static char get_ascii_from_luma(uint8_t luma, const char *scale,
				uint8_t *next_char)
{
	float factor = luma / 256.0;
	uint8_t ascii_id = (uint8_t)(factor * 10);
	char c = scale[*next_char * 10 + ascii_id];
	*next_char = (*next_char + 1) % 9;
	return c;
}
