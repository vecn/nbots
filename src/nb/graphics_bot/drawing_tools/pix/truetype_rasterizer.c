#include <stdio.h>
#include <stdint.h>
#include <string.h>

#define STBTT_STATIC
#define STB_TRUETYPE_IMPLEMENTATION
#include "tiny_libs/stb_truetype.h"

#include "tiny_libs/stb_easy_font.h"

#include "truetype_rasterizer.h"


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
	while (string[c] != '\0') {
		int codepoint = string[c];
		int x0, y0;
		int x1, y1;
		stbtt_GetCodepointBitmapBox(font, codepoint, scale, scale,
					    &x0, &y0, &x1, &y1);
		*w += x1 - x0;
		*h += y1 - y0;
	}
	/* AQUI VOY Soportar multilínea */
}

static void get_size_bitmap_font(const char *ttf_memblock,
				 const char *string, uint16_t size,
				 int *w, int *h)
{
	*w = stb_easy_font_width(string);
	*h = stb_easy_font_height(string);
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

	int x0, y0;
	int x1, y1;
	stbtt_GetFontBoundingBox(font, &x0, &y0, &x1, &y1);

	int bitmap_width = (x1 - x0) * strlen(string);/* AQUI VOY */

	int c = 0;
	while (string[c] != '\0') {
		int codepoint = string[c];
		stbtt_GetCodepointBitmapBox(font, codepoint, scale, scale,
					    &x0, &y0, &x1, &y1);
		int w = x1 - x0;
		int h = y1 - y0;
		stbtt_MakeCodepointBitmap(font, bitmap, w, h, bitmap_width,
					  scale, scale, codepoint);
	}
}

static void bake_bitmap_font(const char *string, uint16_t size,
			     uint8_t *bitmap)
{
	stb_easy_font_print();/* AQUI VOY */
}
