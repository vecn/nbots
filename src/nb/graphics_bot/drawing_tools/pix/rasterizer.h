#ifndef __NB_GRAPHICS_BOT_DRAWING_TOOLS_PIX_RASTERIZER_H__
#define __NB_GRAPHICS_BOT_DRAWING_TOOLS_PIX_RASTERIZER_H__

#include <stdint.h>
#include <stdbool.h>

void nb_graphics_rasterizer_line(int x0, int y0, int x1, int y1,
				 bool antialiasing,
				 void (*set_pixel)(int x, int y, 
						   uint8_t i, void*),
				 void *pixel_data);

void nb_graphics_rasterizer_quad_bezier(int x0, int y0, int x1, int y1,
					int xcontrol, int ycontrol,
					bool antialiasing,
					void (*set_pixel)(int x, int y, 
							  uint8_t i, void*),
					void *pixel_data);

void nb_graphics_rasterizer_quad_rational_bezier(int x0, int y0,
						 int x1, int y1,
						 int xcontrol, int ycontrol,
						 float w,
						 bool antialiasing,
						 void (*set_pixel)
						 	(int x, int y,
							 uint8_t i, void*),
						 void *pixel_data);

void nb_graphics_rasterizer_cubic_bezier(int x0, int y0, int x1, int y1,
					 int x0control, int y0control,
					 int x1control, int y1control,
					 bool antialiasing,
					 void (*set_pixel)(int x, int y, 
							   uint8_t i, void*),
					 void *pixel_data);

#endif
