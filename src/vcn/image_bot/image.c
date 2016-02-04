#include <stdint.h>
#include <string.h>

#include "vcn/image_bot/image.h"

struct vcn_image_s {
	uint32_t width;
	uint32_t height;
	uint8_t comp_x_pixel;
	uint8_t* pixels;
};

inline vcn_image_t* vcn_image_create(void)
{
	return calloc(1, sizeof(vcn_image_t));
}

void vcn_image_destroy(vcn_image_t *img)
{
	if (NULL != img->pixels)
		stbi_image_free(img->pixels);
	free(img);
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
