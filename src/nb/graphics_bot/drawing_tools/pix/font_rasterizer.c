#define STBTT_STATIC
#define STB_TRUETYPE_IMPLEMENTATION

#include "tiny_libs/stb_truetype.h"

#include "font_rasterizer.h"

void fontarara()
{
	char ttf_buffer[1<<25];
/*
	stbtt_fontinfo font;
	int c = (argc > 1 ? atoi(argv[1]) : 'a');
	FILE *fp = fopen("c:/windows/fonts/arialbd.ttf", "rb");
	fread(ttf_buffer, 1, 1<<25, fp);
	int idx = stbtt_GetFontOffsetForIndex(ttf_buffer,0);
	stbtt_InitFont(&font, ttf_buffer, idx);
	float pixels_height = 10;
	float scale = stbtt_ScaleForPixelHeight(&font, pixels_height);
	int w, h;
	unsigned char *bitmap = stbtt_GetCodepointBitmap(&font, 0, scale,
							 c, &w, &h, 0,0);
	for (int j=0; j < h; ++j) {
		for (int i=0; i < w; ++i)
			putchar(" .:ioVM@"[bitmap[j*w+i]>>5]);
		putchar('\n');
	}
*/
}
