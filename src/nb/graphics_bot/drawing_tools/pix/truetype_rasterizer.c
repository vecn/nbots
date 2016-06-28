#include <stdio.h>
#include <stdint.h>
#include <string.h>

#define NB_GRAPHICS_PIX_LETTER_SPACING 0
#define NB_GRAPHICS_PIX_LEADING 125

#define STBTT_STATIC
#define STB_TRUETYPE_IMPLEMENTATION
#include "tiny_libs/stb_truetype.h"

#include "truetype_rasterizer.h"

#define MAX(a,b) (((a)>(b))?(a):(b))

static char* read_ttf_memblock(const char *type);
static const char* get_font_url(const char *type);
static void get_size_truetype_font(const char *ttf_memblock,
				   const char *string, uint16_t size,
				   int *w, int *h);
static void get_size_bitmap_font(const char *ttf_memblock,
				 const char *string, uint16_t size,
				 int *w, int *h);
static void bake_truetype_font(const char *ttf_memblock,
			       const char *string, uint16_t size,
			       uint8_t *bitmap);
static void bake_bitmap_font(const char *string, uint16_t size,
			     uint8_t *bitmap);

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
	stbtt_fontinfo *font = alloca(sizeof(sbtt_fontinfo));
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
		}
		c += 1;
	}
	int spaces = (max_col - 1) * NB_GRAPHICS_PIX_LETTER_SPACING;
	*w = max_w + spaces;
	float leading = NB_GRAPHICS_PIX_LEADING / 100.0f;
	int line_spaces = (N_rows - 1) * (int)(leading * size);
	*h = size * N_rows + line_spaces;
}

static void get_size_bitmap_font(const char *ttf_memblock,
				 const char *string, uint16_t size,
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
		}
		col += 1;
		c += 1;
	}
	int spaces = (max_col - 1) * NB_GRAPHICS_PIX_LETTER_SPACING;
	*w = size * max_col + spaces;
	float leading = NB_GRAPHICS_PIX_LEADING / 100.0f;
	int line_spaces = (N_rows - 1) * (int)(leading * size);
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
	FILE *fp = fopen(get_font_url(type), "rb");
	char *ttf_memblock = NULL;
	if (NULL != fp) {
		fseek(fp, 0L, SEEK_END);
		int64_t file_size  = ftell(fp);
		fseek(fp, 0L, SEEK_SET);

		ttf_memblock = malloc(file_size);
		fread(ttf_memblock, 1, ttf_memblock, fp);
		fclose(fp);
	}
	return ttf_memblock;
}

static const char* get_font_url(const char *type)
{
	return "/home/victor/repos/nbots/src/nb/graphics_bot/drawing_tools/pix/fonts/FreeSans.ttf";
}

static void bake_truetype_font(const char *ttf_memblock,
			       const char *string, uint16_t size,
			       uint8_t *bitmap)
{
	stbtt_fontinfo *font = alloca(sizeof(sbtt_fontinfo));
	int font_idx = stbtt_GetFontOffsetForIndex(ttf_memblock,0);
	stbtt_InitFont(font, ttf_memblock, font_idx);
	float scale = stbtt_ScaleForPixelHeight(font, size);

	int bm_w, bm_h;
	get_size_truetype_font(ttf_memblock, string, size, &bm_w, &bm_h);

	int c = 0;
	while ('\0' != string[c]) {
		int codepoint = string[c];
		int x0, y0;
		int x1, y1;
		stbtt_GetCodepointBitmapBox(font, codepoint, scale, scale,
					    &x0, &y0, &x1, &y1);
		int w = x1 - x0;
		int h = y1 - y0;
		stbtt_MakeCodepointBitmap(font, bitmap, w, h, bitmap_width,
					  scale, scale, codepoint);
		c += 1;
	}
	/* AQUI VOY */
}

static void bake_bitmap_font(const char *string, uint16_t size,
			     uint8_t *bitmap)
{
	int bm_w, bm_h;
	get_size_bitmap_font(string, size, &bm_w, &bm_h);

	char **font8x8 = font8x8_basic; 

	int col = 0;
	int row = 0;
	int c = 0;
	while ('\0' != string[c]) {
		if ('\n' == string[c]) {
			col = 0;
			row += 1;
		} else {
			col += 1;
			int codepoint = string[c] % 128;
			char char8x8[8] = font8x8[codepoint];
/* AQUI VOY */
			if (8 > size) {
			} else if (8 == size) {
				
			} else {

			}
		}
		c += 1;
	}
}
