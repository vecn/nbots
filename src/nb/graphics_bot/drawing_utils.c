#include "nb/graphics_bot/drawing_utils.h"

void nb_drawing_utils_set_center_and_zoom(camera_t *cam, const double box[4],
					  double width, double height)
{
	cam->center[0] = (box[0] + box[2]) / 2.0;
	cam->center[1] = (box[1] + box[3]) / 2.0;
	cam->zoom = width / (box[2] - box[0]);
	if (cam->zoom > height / (box[3] - box[1]))
		cam->zoom = height / (box[3] - box[1]);
	cam->zoom *= 0.9;
	cam->width = width;
	cam->height = height;
}
