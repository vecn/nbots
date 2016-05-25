#ifndef __NB_GEOMETRIC_BOT_MESH_MODULES2D_EXPORTER_CAIRO_DRAWING_UTILS_H__
#define __NB_GEOMETRIC_BOT_MESH_MODULES2D_EXPORTER_CAIRO_DRAWING_UTILS_H__

static double _NB_RGB_BLACK[4] = {0.0, 0.0, 0.0};
static double _NB_RGB_RED[3] = {0.0, 0.0, 1.0};
static double _NB_RGB_BLUE[3] = {0.0, 0.0, 1.0};

typedef struct {
	int width;
	int height;
	double center[2];
	double zoom;
} camera_t;

void nb_drawing_utils_set_center_and_zoom(camera_t *cam, const double box[4],
					  double width, double height);

#endif
