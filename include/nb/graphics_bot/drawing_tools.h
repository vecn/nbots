#ifndef __NB_GRAPHICS_BOT_DRAWING_TOOLS_H__
#define __NB_GRAPHICS_BOT_DRAWING_TOOLS_H__

#include <stdbool.h>
#include <stdint.h>

typedef struct nb_graphics_context_s nb_graphics_context_t;
typedef struct nb_graphics_palette_s nb_graphics_palette_t;

typedef enum {
	NB_LINEAR, NB_RADIAL
} nb_graphics_grad_t;

typedef enum {
	NB_RAINBOW,
	NB_SUNSET,
	NB_FRENCH
} nb_graphics_palette_preset;

typedef enum {
	NB_WHITE, NB_BLACK,
	NB_GRAY, NB_LIGHT_GRAY, NB_DARK_GRAY,
	NB_RED, NB_GREEN, NB_BLUE,
	NB_YELLOW, NB_CYAN, NB_MAGENTA,
	NB_ORANGE, NB_AQUAMARIN, NB_VIOLET,
	NB_ROSE, NB_CHARTREUSE, NB_AZURE
} nb_graphics_color_t;

typedef struct {
	float width;
	float height;
	float x_left;
	float y_top;
} nb_graphics_text_attr_t;

typedef struct {
	int width;
	int height;
	float center[2];
	float zoom;
} nb_graphics_camera_t;

void nb_graphics_export(const char* filename, int width, int height,
			void (*draw)(nb_graphics_context_t *g, int w, int h,
				     const void *const data),
			const void *const data);

nb_graphics_camera_t* nb_graphics_get_camera(nb_graphics_context_t *g);
bool nb_graphics_is_camera_enabled(const nb_graphics_context_t *g);
void nb_graphics_disable_camera(nb_graphics_context_t *g);
void nb_graphics_enable_camera(nb_graphics_context_t *g);

void nb_graphics_move_to(nb_graphics_context_t *g, float x, float y);
void nb_graphics_line_to(nb_graphics_context_t *g, float x, float y);
void nb_graphics_qcurve_to(nb_graphics_context_t *g, float x, float y,
			   float xcontrol, float ycontrol);
void nb_graphics_qrcurve_to(nb_graphics_context_t *g, float x, float y,
			    float xcontrol, float ycontrol, float w);
void nb_graphics_curve_to(nb_graphics_context_t *g, float x, float y,
			  float x0_control, float y0_control,
			  float x1_control, float y1_control);

void nb_graphics_set_ellipse(nb_graphics_context_t *g, float x, float y,
			     float rx, float ry, float angle);

void nb_graphics_set_circle(nb_graphics_context_t *g,
			    float x, float y, float r);

void nb_graphics_set_point(nb_graphics_context_t *g,
			   float x, float y, float size);

void nb_graphics_set_rectangle(nb_graphics_context_t *g, float x1, float y1,
				  float x2, float y2);

void nb_graphics_close_path(nb_graphics_context_t *g);

void nb_graphics_set_line_width(nb_graphics_context_t *g, float w);

void nb_graphics_set_source(nb_graphics_context_t *g,
			    nb_graphics_color_t color);
void nb_graphics_set_source_rgb(nb_graphics_context_t *gctx,
				uint8_t r, uint8_t g, uint8_t b);

void nb_graphics_set_source_rgba(nb_graphics_context_t *gctx,
				 uint8_t r, uint8_t g, uint8_t b, uint8_t a);

void nb_graphics_set_source_grad(nb_graphics_context_t *g,
				 nb_graphics_grad_t grad,
				 float x1, float y1,
				 float x2, float y2,
				 nb_graphics_palette_t *pat);

void nb_graphics_set_source_trg(nb_graphics_context_t *g,
				float x1, float y1,
				float x2, float y2,
				float x3, float y3,
				const uint8_t rgba1[4],
				const uint8_t rgba2[4],
				const uint8_t rgba3[4]);

void nb_graphics_fill(nb_graphics_context_t *g);

void nb_graphics_fill_preserve(nb_graphics_context_t *g);

void nb_graphics_stroke(nb_graphics_context_t *g);

void nb_graphics_set_font_type(nb_graphics_context_t *g, const char *type);

void nb_graphics_set_font_size(nb_graphics_context_t *g, uint16_t size);

void nb_graphics_show_text(nb_graphics_context_t *g, const char *str);
void nb_graphics_get_text_attr(const nb_graphics_context_t *g, const char *label,
			       nb_graphics_text_attr_t *attr);

void nb_graphics_cam_fit_box(nb_graphics_camera_t *cam, const double box[4],
			     float width, float height);

nb_graphics_palette_t* nb_graphics_palette_create(void);
nb_graphics_palette_t* nb_graphics_palette_create_preset(nb_graphics_palette_preset preset);
void nb_graphics_palette_destroy(nb_graphics_palette_t* palette);
void nb_graphics_palette_clear(nb_graphics_palette_t* palette);
void nb_graphics_palette_add_rgba(nb_graphics_palette_t* palette, float tic,
				  uint8_t r, uint8_t g, uint8_t b, uint8_t a);
void nb_graphics_palette_get_rgba(const nb_graphics_palette_t *const palette,
				  float factor,
				  uint8_t rgba[4]);

void nb_graphics_palette_draw(nb_graphics_context_t *g, 
			      const nb_graphics_palette_t *const palette,
			      float x, float y, float w, float h,
			      float border, float min_v, float max_v);

#endif
