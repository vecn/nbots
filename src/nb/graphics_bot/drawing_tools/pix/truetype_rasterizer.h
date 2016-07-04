#ifndef __NB_GRAPHICS_BOT_DRAWING_TOOLS_PIX_TRUETYPE_RASTERIZER_H__
#define __NB_GRAPHICS_BOT_DRAWING_TOOLS_PIX_TRUETYPE_RASTERIZER_H__

void nb_graphics_truetype_rasterizer_get_size(const char *string,
					      const char *type, uint16_t size,
					      int *w, int *h);
void nb_graphics_truetype_rasterizer_bake(const char *string,
					  const char *type, uint16_t size,
					  uint8_t *bitmap);

#endif
