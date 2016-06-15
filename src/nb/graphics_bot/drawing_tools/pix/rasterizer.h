#ifndef __NB_GRAPHICS_BOT_DRAWING_TOOLS_PIX_RASTERIZER_H__
#define __NB_GRAPHICS_BOT_DRAWING_TOOLS_PIX_RASTERIZER_H__


void bresenham_line(int x0, int y0, int x1, int y1,
		    void (*set_pixel)(void*),
		    void *set_pixel_data);

void bresenham_ellipse(int xc, int yc, int rx, int ry,
		       void (*set_pixel)(void*),
		       void *set_pixel_data);

void bresenham_circle(int xc, int yc, int r,
		      void (*set_pixel)(void*),
		      void *set_pixel_data);

#endif
