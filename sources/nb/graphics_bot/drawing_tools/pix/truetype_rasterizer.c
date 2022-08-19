#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "nb/memory_bot.h"

#define NB_GRAPHICS_PIX_LETTER_SPACING 0
#define NB_GRAPHICS_PIX_LEADING 125

#ifndef NB_EXCLUDE_IMPORTED_LIBS
    #define STB_TRUETYPE_IMPLEMENTATION
    #include "imported_libs/stb_truetype.h"
#else
    typedef int stbtt_fontinfo;
#endif

#include "font_bitmap.h"

#include "truetype_rasterizer.h"

#define MAX(a,b) (((a)>(b))?(a):(b))


#define GET_FONT_URL(file, type)					\
	do {								\
		int type_len = strlen(type);				\
		char *path = getenv("NB_FONTS_DIR");			\
		if (NULL == path) {					\
			file = nb_allocate_on_stack(type_len + 10);			\
			sprintf(file, "/home/.fonts/%s.ttf", type);	\
		} else {						\
			int path_len = strlen(path);			\
			if ('/' == path[path_len - 1]) {		\
				file = nb_allocate_on_stack(path_len + type_len + 4);	\
				sprintf(file, "%s%s.ttf", path, type);	\
			} else {					\
				file = nb_allocate_on_stack(path_len + type_len + 5);	\
				sprintf(file, "%s/%s.ttf", path, type);	\
			}						\
		}							\
	} while(0)

static unsigned char* read_ttf_memblock(const char *type);
static void get_size_truetype_font(stbtt_fontinfo *font,
				   const char *string, uint16_t size,
				   int *w, int *h);
static void get_size_bitmap_font(const char *string, uint16_t size,
				 int *w, int *h);
static void bake_truetype_font(const unsigned char *ttf_memblock,
			       const char *string, uint16_t size,
			       uint8_t *bitmap);
static void bake_bitmap_font(const char *string, uint16_t size,
			     uint8_t *bitmap);
static void stamp_bmchar_to_bitmap(uint8_t *bitmap, uint16_t size,
				   const char char8x8[8], int col,
				   int row, int bm_width);
static uint8_t bmchar_get_intensity(const char char8x8[8], uint16_t size,
				    int i, int j);
static float bmchar_get_intensity_eq8(const char char8x8[8],
				      int i, int j);
static float bmchar_get_intensity_lt8(const char char8x8[8], uint16_t size,
				      int i, int j);
static float bmchar_get_intensity_gt8(const char char8x8[8], uint16_t size,
				      int i, int j);

void nb_graphics_truetype_rasterizer_get_size(const char *string,
					      const char *type, uint16_t size,
					      int *w, int *h)
{
#ifndef NB_EXCLUDE_IMPORTED_LIBS
	unsigned char *ttf_memblock = read_ttf_memblock(type);
	if (NULL != ttf_memblock) {
		stbtt_fontinfo *font = nb_allocate_on_stack(sizeof(stbtt_fontinfo));
		int font_idx = stbtt_GetFontOffsetForIndex(ttf_memblock, 0);
		stbtt_InitFont(font, ttf_memblock, font_idx);
		get_size_truetype_font(font, string, size, w, h);
		nb_free_mem(ttf_memblock);
	} else {
		get_size_bitmap_font(string, size, w, h);
	}
#else
	*w = size * strlen(string);
	*h = size;
#endif
}

static void get_size_truetype_font(stbtt_fontinfo *font,
				   const char *string, uint16_t size,
				   int *w, int *h)
{
#ifndef NB_EXCLUDE_IMPORTED_LIBS
	float scale = stbtt_ScaleForPixelHeight(font, size);

	*w = 0;
	*h = 0;
	int c = 0;
	int N_rows = 1;
	int max_w = 0;
	while ('\0' != string[c]) {
		int letter_spacing = NB_GRAPHICS_PIX_LETTER_SPACING;
		int codepoint = string[c];
		int advance, lsb;
		stbtt_GetCodepointHMetrics(font, codepoint,
					   &advance, &lsb);
		if ('\n' == string[c]) {
			N_rows += 1;
			max_w = MAX(max_w, *w);
			*w = 0;
		} else {
			int x0, y0;
			int x1, y1;
			stbtt_GetCodepointBitmapBox(font, codepoint,
						    scale, scale,
						    &x0, &y0, &x1, &y1);
			*w += scale * advance + letter_spacing;
		}
		c += 1;
	}
	max_w = MAX(max_w, *w);

	*w = max_w;
	float leading = NB_GRAPHICS_PIX_LEADING / 100.0f;
	int line_spaces = N_rows * ((int)(leading * size) - size);
	*h = size * N_rows + line_spaces;
#endif
}

static void get_size_bitmap_font(const char *string, uint16_t size,
				 int *w, int *h)
{
	int max_col = 0;
	int N_rows = 1;
	int col = 0;
	int c = 0;
	while ('\0' != string[c]) {
		if ('\n' == string[c]) {
			N_rows += 1;
			max_col = MAX(max_col, col);
			col = 0;
		} else {
			col += 1;
		}
		c += 1;
	}
	max_col = MAX(max_col, col);

	int spaces = (max_col - 1) * NB_GRAPHICS_PIX_LETTER_SPACING;
	*w = size * max_col + spaces;
	float leading = NB_GRAPHICS_PIX_LEADING / 100.0f;
	int line_spaces = (N_rows - 1) * ((int)(leading * size) - size);
	*h = size * N_rows + line_spaces;
}

void nb_graphics_truetype_rasterizer_bake(const char *string,
					  const char *type, uint16_t size,
					  uint8_t *bitmap)
{
	unsigned char *ttf_memblock = read_ttf_memblock(type);
	if (NULL != ttf_memblock) {
		bake_truetype_font(ttf_memblock, string, size, bitmap);
		nb_free_mem(ttf_memblock);
	} else {
		bake_bitmap_font(string, size, bitmap);
	}
}

static unsigned char* read_ttf_memblock(const char *type)
{
	char *filename;
	GET_FONT_URL(filename, type);
	FILE *fp = fopen(filename, "rb");
	unsigned char *ttf_memblock = NULL;
	if (NULL != fp) {
		fseek(fp, 0L, SEEK_END);
		int64_t file_size  = ftell(fp);
		fseek(fp, 0L, SEEK_SET);

		ttf_memblock = nb_allocate_mem(file_size);
		fread(ttf_memblock, 1, file_size, fp);
		fclose(fp);
	}
	return ttf_memblock;
}

static void bake_truetype_font(const unsigned char *ttf_memblock,
			       const char *string, uint16_t size,
			       uint8_t *bitmap)
{
	stbtt_fontinfo *font = nb_allocate_on_stack(sizeof(stbtt_fontinfo));

#ifndef NB_EXCLUDE_IMPORTED_LIBS
	int font_idx = stbtt_GetFontOffsetForIndex(ttf_memblock,0);
	stbtt_InitFont(font, ttf_memblock, font_idx);
	float scale = stbtt_ScaleForPixelHeight(font, size);
#endif

      	int bm_w, bm_h;
	get_size_truetype_font(font, string, size, &bm_w, &bm_h);

#ifndef NB_EXCLUDE_IMPORTED_LIBS
	int ascent;
	stbtt_GetFontVMetrics(font, &ascent, 0, 0);
	int baseline = (int) (scale * ascent);

	int line_spacing = (int)(size * (NB_GRAPHICS_PIX_LEADING/100.f));

	int cumulative_width = 0;
	int row = 0;
	int c = 0;
	while ('\0' != string[c]) {
		int letter_spacing = NB_GRAPHICS_PIX_LETTER_SPACING;
		int codepoint = string[c];
		int advance, lsb;
		stbtt_GetCodepointHMetrics(font, codepoint,
					   &advance, &lsb);
		if ('\n' == string[c]) {
			cumulative_width = 0;
			row += 1;
		} else {
			int x0, y0;
			int x1, y1;
			stbtt_GetCodepointBitmapBox(font, codepoint,
						    scale, scale,
						    &x0, &y0, &x1, &y1);
			int w = x1 - x0;
			int h = y1 - y0;
			int mean_offset = MAX(0, baseline + y0);
			int pix_row = row * line_spacing + mean_offset;
			int pix_col = cumulative_width;
			stbtt_MakeCodepointBitmap(font, 
						  &(bitmap[pix_row * bm_w +
							   pix_col]),
						  w, h, bm_w,
						  scale, scale,
						  codepoint);
			
			cumulative_width += scale * advance + letter_spacing;
		}
		c += 1;
	}
#endif
}

static void bake_bitmap_font(const char *string, uint16_t size,
			     uint8_t *bitmap)
{
	int bm_w, bm_h;
	get_size_bitmap_font(string, size, &bm_w, &bm_h);

	char *font8x8 = (char*) nb_font8x8_basic; 

	int col = 0;
	int row = 0;
	int c = 0;
	while ('\0' != string[c]) {
		if ('\n' == string[c]) {
			col = 0;
			row += 1;
		} else {
			int codepoint = string[c] % 128;
			char char8x8[8];
			memcpy(char8x8, font8x8 + 8 * codepoint, 8);
			stamp_bmchar_to_bitmap(bitmap, size, char8x8,
					       col, row, bm_w);
			col += 1;
		}
		c += 1;
	}
}

static void stamp_bmchar_to_bitmap(uint8_t *bitmap, uint16_t size,
				   const char char8x8[8], int col,
				   int row, int bm_width)
{
	int line_spacing = (int)(size * (NB_GRAPHICS_PIX_LEADING/100.f));
	int letter_spacing = NB_GRAPHICS_PIX_LETTER_SPACING;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			int pix_row = row * line_spacing + i;
			int pix_col = col * (size + letter_spacing) + j;
			bitmap[pix_row * bm_width + pix_col] = 
				bmchar_get_intensity(char8x8, size, i, j);
		}
	}
}

static uint8_t bmchar_get_intensity(const char char8x8[8], uint16_t size,
				    int i, int j)
{
	float intensity;
	if (8 == size)
		intensity = bmchar_get_intensity_eq8(char8x8, i, j);
	else if (8 > size)
		intensity = bmchar_get_intensity_lt8(char8x8, size, i, j);
	else
		intensity = bmchar_get_intensity_gt8(char8x8, size, i, j);
	return (int)(intensity * 255 + 0.5);
}

static float bmchar_get_intensity_eq8(const char char8x8[8],
				      int i, int j)
{
	char bit_row = char8x8[i];
	char bit_mask = 1 << j;
	bool bit_enabled = (0 != (bit_row & bit_mask));
	return (bit_enabled)?(1.0f):(0.0f);
}

static float bmchar_get_intensity_lt8(const char char8x8[8], uint16_t size,
				      int i, int j)
{
	float pixels = 8.0f / size;
	int ifirst = (int)(i * pixels);
	int ilast = (int)((i+1) * pixels);
	if ((i+1) * pixels - ilast > 0.0)
		ilast += 1;
	int jfirst = (int)(j * pixels);
	int jlast = (int)((j+1) * pixels);
	if ((j+1) * pixels - jlast > 0.0)
		jlast += 1;

	float weight = 0.0f;
	float total = 0;
	for (int k = ifirst; k < ilast; k++) {
		char bit_row = char8x8[k];
		float yw = 1.0f;
		if (k == ifirst) {
			float rb = ifirst + 1;
			yw = (rb - i * pixels);
		} else if (k == ilast - 1) {
			float rb = ilast;
			yw = (1.0f - rb + (i+1) * pixels);
		}
		for (int l = jfirst; l < jlast; l++) {
			char bit_mask = 1 << l;
			bool bit_enabled = (0 != (bit_row & bit_mask));
				
			float xw = 1.0f;
			if (l == jfirst) {
				float rb = jfirst + 1;
				xw = (rb - j * pixels);
			} else if (l == jlast - 1) {
				float rb = jlast;
				xw = (1.0f - rb + (j+1) * pixels);
			}

			if (bit_enabled)
				weight += xw * yw;
				
			total += xw * yw;
		}
	}
	weight /= total;
	return weight;
}

static float bmchar_get_intensity_gt8(const char char8x8[8], uint16_t size,
				      int i, int j)
{
	float pixels = size / 8.0f;

	float weight = 1.0f;

	int k = (int)(i/pixels);
	if (i/pixels - k + 1.0f/pixels >= 1.0)
		weight *= i/pixels - k;

	int l = (int)(j/pixels);
	if (j/pixels - l + 1.0/pixels >= 1.0)
		weight *= j/pixels - l;
	
	char bit_row = char8x8[k];
	char bit_mask = 1 << l;
	bool bit_enabled = (0 != (bit_row & bit_mask));
	
	return (bit_enabled)?(1.0f):(0.0f);
}
