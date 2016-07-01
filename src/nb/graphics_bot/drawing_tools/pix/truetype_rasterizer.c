#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <alloca.h>
#include <stdbool.h>

#define NB_GRAPHICS_PIX_LETTER_SPACING 0
#define NB_GRAPHICS_PIX_LEADING 125

#define STBTT_STATIC
#define STB_TRUETYPE_IMPLEMENTATION
#include "tiny_libs/stb_truetype.h"

#include "font_bitmap.h"

#include "truetype_rasterizer.h"

#define MAX(a,b) (((a)>(b))?(a):(b))


#define GET_FONT_URL(file, type)					\
	do {								\
		int type_len = strlen(type);				\
		char *path = getenv("NB_FONTS_DIR");			\
		if (NULL == path) {					\
			file = alloca(type_len + 10);			\
			sprintf(file, "/home/.fonts/%s.ttf", type);	\
		} else {						\
			int path_len = strlen(path);			\
			if ('/' == path[path_len - 1]) {		\
				file = alloca(path_len + type_len + 4);	\
				sprintf(file, "%s%s.ttf", path, type);	\
			} else {					\
				file = alloca(path_len + type_len + 5);	\
				sprintf(file, "%s/%s.ttf", path, type);	\
			}						\
		}							\
	} while(0)

#define INIT_CHAR_BITMAP(char_bitmap)					\
	do {								\
		if (100 < size)						\
			char_bitmap = malloc(size * size);		\
		else							\
			char_bitmap = alloca(size * size);		\
		memset(char_bitmap, 0, size * size);			\
	} while(0)

#define FINISH_CHAR_BITMAP(char_bitmap)					\
	do {								\
		if (100 < size)						\
			free(char_bitmap);				\
	} while(0)

static char* read_ttf_memblock(const char *type);
static void get_size_truetype_font(const char *ttf_memblock,
				   const char *string, uint16_t size,
				   int *w, int *h);
static void get_size_bitmap_font(const char *string, uint16_t size,
				 int *w, int *h);
static void bake_truetype_font(const char *ttf_memblock,
			       const char *string, uint16_t size,
			       uint8_t *bitmap);
static void stamp_tt_bmchar_to_bitmap(uint8_t *bitmap, int w, int h,
				      char *char_bitmap, int col,
				      int row, int bm_width);
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
	char *ttf_memblock = read_ttf_memblock(type);
	if (NULL != ttf_memblock) {
		get_size_truetype_font(ttf_memblock, string, size, w, h);
		free(ttf_memblock);
	} else {
		get_size_bitmap_font(string, size, w, h);
	}
}

static void get_size_truetype_font(const char *ttf_memblock,
				   const char *string, uint16_t size,
				   int *w, int *h)
{
	stbtt_fontinfo *font = alloca(sizeof(stbtt_fontinfo));
	int font_idx = stbtt_GetFontOffsetForIndex(ttf_memblock,0);
	stbtt_InitFont(font, ttf_memblock, font_idx);
	float scale = stbtt_ScaleForPixelHeight(font, size);

	*w = 0;
	*h = 0;
	int c = 0;
	int N_rows = 1;
	int max_w = 0;
	int max_col = 0;
	int col = 0;
	while ('\0' != string[c]) {
		if ('\n' == string[c]) {
			N_rows += 1;
			max_w = MAX(max_w, *w);
			*w = 0;
			max_col = MAX(max_col, col);
			col = 0;
		} else {
			int codepoint = string[c];
			int x0, y0;
			int x1, y1;
			stbtt_GetCodepointBitmapBox(font, codepoint,
						    scale, scale,
						    &x0, &y0, &x1, &y1);
			*w += x1 - x0;
			col += 1;
			printf("HEIGHT: %i\n", y1-y0);/**/
		}
		c += 1;
	}
	max_w = MAX(max_w, *w);
	max_col = MAX(max_col, col);

	int spaces = (max_col - 1) * NB_GRAPHICS_PIX_LETTER_SPACING;
	*w = max_w + spaces;
	float leading = NB_GRAPHICS_PIX_LEADING / 100.0f;
	int line_spaces = (N_rows - 1) * ((int)(leading * size) - size);
	*h = size * N_rows + line_spaces;
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
	char *ttf_memblock = read_ttf_memblock(type);
	if (NULL != ttf_memblock) {
		bake_truetype_font(ttf_memblock, string, size, bitmap);
		free(ttf_memblock);
	} else {
		bake_bitmap_font(string, size, bitmap);
	}
}

static char* read_ttf_memblock(const char *type)
{
	char *filename;
	GET_FONT_URL(filename, type);
	FILE *fp = fopen(filename, "rb");
	char *ttf_memblock = NULL;
	if (NULL != fp) {
		fseek(fp, 0L, SEEK_END);
		int64_t file_size  = ftell(fp);
		fseek(fp, 0L, SEEK_SET);

		ttf_memblock = malloc(file_size);
		fread(ttf_memblock, 1, file_size, fp);
		fclose(fp);
	}
	return ttf_memblock;
}

static void bake_truetype_font(const char *ttf_memblock,
			       const char *string, uint16_t size,
			       uint8_t *bitmap)
{
	stbtt_fontinfo *font = alloca(sizeof(stbtt_fontinfo));
	int font_idx = stbtt_GetFontOffsetForIndex(ttf_memblock,0);
	stbtt_InitFont(font, ttf_memblock, font_idx);
	float scale = stbtt_ScaleForPixelHeight(font, size);

	int bm_w, bm_h;
	get_size_truetype_font(ttf_memblock, string, size, &bm_w, &bm_h);

	char* char_bitmap;
	INIT_CHAR_BITMAP(char_bitmap);

	int col = 0;
	int row = 0;
	int c = 0;
	while ('\0' != string[c]) {
		if ('\n' == string[c]) {
			col = 0;
			row += 1;
		} else {
			int codepoint = string[c];
			int x0, y0;
			int x1, y1;
			stbtt_GetCodepointBitmapBox(font, codepoint,
						    scale, scale,
						    &x0, &y0, &x1, &y1);
			int w = x1 - x0;
			int h = y1 - y0;
			stbtt_MakeCodepointBitmap(font, bitmap, w, h,
						  w, scale, scale,
						  codepoint);
			stamp_tt_bmchar_to_bitmap(bitmap, w, h, char_bitmap,
						  col, row, bm_w);
			col += 1;
		}
		c += 1;
	}/*
	while (text[ch]) {
		int advance,lsb,x0,y0,x1,y1;
		float x_shift = xpos - (float) floor(xpos);
		stbtt_GetCodepointHMetrics(&font, text[ch], &advance, &lsb);
		stbtt_GetCodepointBitmapBoxSubpixel(&font, text[ch], scale,scale,x_shift,0, &x0,&y0,&x1,&y1);
		stbtt_MakeCodepointBitmapSubpixel(&font, &screen[baseline + y0][(int) xpos + x0], x1-x0,y1-y0, 79, scale,scale,x_shift,0, text[ch]);
		xpos += (advance * scale);
		if (text[ch+1])
			xpos += scale*stbtt_GetCodepointKernAdvance(&font, text[ch],text[ch+1]);
		++ch;
	}*/
	FINISH_CHAR_BITMAP(char_bitmap);
}

static void stamp_tt_bmchar_to_bitmap(uint8_t *bitmap, int w, int h,
				      char *char_bitmap, int col,
				      int row, int bm_width)
{
	int line_spacing = (int)(h * (NB_GRAPHICS_PIX_LEADING/100.f));
	int letter_spacing = NB_GRAPHICS_PIX_LETTER_SPACING;
	for (int i = 0; i < h; i++) {
		for (int j = 0; j < w; j++) {
			int pix_row = row * line_spacing + i;
			int pix_col = col * (w + letter_spacing) + j;
			bitmap[pix_row * bm_width + pix_col] = 
				char_bitmap[i * w + j];
		}
	}	
}

static void bake_bitmap_font(const char *string, uint16_t size,
			     uint8_t *bitmap)
{
	int bm_w, bm_h;
	get_size_bitmap_font(string, size, &bm_w, &bm_h);

	char *font8x8 = font8x8_basic; 

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
