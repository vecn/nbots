/* Rasterizer based on Bresenham algorithm */

#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include "rasterizer.h"

#define PIX_STACK_STATIC_SIZE 500
#define PIX_STACK_DYNAMIC_INCREMENTS 500

#define POW2(a) ((a)*(a))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

typedef struct {
	uint32_t w;
	uint32_t N;
	uint32_t N_alloc;
	uint32_t static_mem[PIX_STACK_STATIC_SIZE];
	uint32_t *dynamic_mem;
} pix_stack_t;

static void line(int x0, int y0, int x1, int y1,
		 void (*set_pixel)(int x, int y,
				   uint8_t i, void*),
		 void *pixel_data);
static void quad_bezier_sgm(int x0, int y0, int x1, int y1,
			    int xcontrol, int ycontrol,
			    void (*set_pixel)(int x, int y,
					      uint8_t i, void*),
			    void *pixel_data);
static void quad_bezier(int x0, int y0, int x1, int y1,
			int xcontrol, int ycontrol, 
			void (*set_pixel)(int x, int y,
					  uint8_t i, void*),
			void *pixel_data,
			void (*qb_sgm)(int, int, int, int, int, int,
				       void (*set_pixel)
				       (int, int, uint8_t, void*),
				       void*));
static void quad_rational_bezier_sgm
				(int x0, int y0, int x1, int y1,
				 int xcontrol, int ycontrol, float w, 
				 void (*set_pixel)(int x, int y,
						   uint8_t i, void*),
				 void *pixel_data);
static void quad_rational_bezier(int x0, int y0, int x1, int y1,
				 int xcontrol, int ycontrol, float w,
				 void (*set_pixel)(int x, int y,
						   uint8_t i, void*),
				 void *pixel_data,
				 void (*qrb_sgm)(int, int, int, int,
						 int, int, float, 
						 void (*set_pixel)
						 (int, int, uint8_t, void*),
						 void *));
static void cubic_bezier_sgm(int x0, int y0, int x1, int y1,
			     float x0_control, float y0_control,
			     float x1_control, float y1_control,
			     void (*set_pixel)(int x, int y,
					       uint8_t i, void*),
			     void *pixel_data);
static void cubic_bezier(int x0, int y0, int x1, int y1,
			 int x0_control, int y0_control,
			 int x1_control, int y1_control,
			 void (*set_pixel)(int x, int y,
					   uint8_t i, void*),
			 void *pixel_data,
			 void (*cb_sgm)(int, int, int, int,
					float, float, float, float,
					void (*set_pixel)
					(int, int, uint8_t, void*),
					void*));
static void line_aa(int x0, int y0, int x1, int y1,
		    void (*set_pixel)(int x, int y,
				      uint8_t i, void*),
		    void *pixel_data);
static void quad_bezier_sgm_aa(int x0, int y0, int x1, int y1,
			       int xcontrol, int ycontrol,
			       void (*set_pixel)(int x, int y,
						 uint8_t i, void*),
			       void *pixel_data);
static void quad_rational_bezier_sgm_aa
				(int x0, int y0, int x1, int y1,
				 int xcontrol, int ycontrol, float w, 
				 void (*set_pixel)(int x, int y,
						   uint8_t i, void*),
				 void *pixel_data);

static void cubic_bezier_sgm_aa(int x0, int y0, int x1, int y1,
				float x0_control, float y0_control,
				float x1_control, float y1_control,
				void (*set_pixel)(int x, int y,
						  uint8_t i, void*),
				void *pixel_data);
static void init_stack(pix_stack_t *ps, int32_t img_width);
static bool pop(pix_stack_t *ps, int *x, int *y);
static void push(pix_stack_t *ps, int x, int y);
static void finish_stack(pix_stack_t *ps);

void nb_graphics_rasterizer_line(int x0, int y0, int x1, int y1,
				 bool antialiasing,
				 void (*set_pixel)(int x, int y, 
						   uint8_t i, void*),
				 void *pixel_data)
{
	if (antialiasing)
		line_aa(x0, y0, x1, y1, set_pixel, pixel_data);
	else
		line(x0, y0, x1, y1, set_pixel, pixel_data);
}

void nb_graphics_rasterizer_quad_bezier(int x0, int y0, int x1, int y1,
					int xcontrol, int ycontrol,
					bool antialiasing,
					void (*set_pixel)(int x, int y, 
							  uint8_t i, void*),
					void *pixel_data)
{
	void (*qb_sgm)(int, int, int, int, int, int,
		       void (*set_pixel)(int, int, uint8_t, void*),
		       void*);
	if (antialiasing)
		qb_sgm = quad_bezier_sgm_aa;
	else
		qb_sgm = quad_bezier_sgm;
	quad_bezier(x0, y0, x1, y1, xcontrol, ycontrol,
		    set_pixel, pixel_data, qb_sgm);

}

void nb_graphics_rasterizer_quad_rational_bezier(int x0, int y0,
						 int x1, int y1,
						 int xcontrol, int ycontrol,
						 float w,
						 bool antialiasing,
						 void (*set_pixel)
						 	(int x, int y,
							 uint8_t i, void*),
						 void *pixel_data)
{
	void (*qrb_sgm)(int, int, int, int, int, int, float, 
			void (*set_pixel)(int, int, uint8_t, void*), 	
			void *);
	if (antialiasing)
		qrb_sgm = quad_rational_bezier_sgm_aa;
	else
		qrb_sgm = quad_rational_bezier_sgm;

	quad_rational_bezier(x0, y0, x1, y1,
			     xcontrol, ycontrol, w,
			     set_pixel, pixel_data,
			     qrb_sgm);
}

void nb_graphics_rasterizer_cubic_bezier(int x0, int y0, int x1, int y1,
					 int x0control, int y0control,
					 int x1control, int y1control,
					 bool antialiasing,
					 void (*set_pixel)(int x, int y, 
							   uint8_t i, void*),
					 void *pixel_data)
{
	void (*cb_sgm)(int, int, int, int, float, float, float, float,
		       void (*set_pixel)(int, int, uint8_t, void*),
		       void*);
	if (antialiasing)
		cb_sgm = cubic_bezier_sgm_aa;
	else
		cb_sgm = cubic_bezier_sgm;

	cubic_bezier(x0, y0, x1, y1,
		     x0control, y0control,
		     x1control, y1control,
		     set_pixel, pixel_data,
		     cb_sgm);
}

static void line(int x0, int y0, int x1, int y1,
		 void (*set_pixel)(int x, int y,
				   uint8_t i, void*),
		 void *pixel_data)
{
	int sx = (x0 < x1)? 1 : -1;
	int sy = (y0 < y1)? 1 : -1;
	int dx =  abs(x1 - x0);
	int dy = -abs(y1 - y0);
	int err = dx + dy;
                                                    
	while (true) {
		set_pixel(x0, y0, 255, pixel_data);
		int e2 = 2 * err;
		if (e2 >= dy) {
			if (x0 == x1)
				break;
			err += dy;
			x0 += sx;
		}
		if (e2 <= dx) {
			if (y0 == y1)
				break;
			err += dx;
			y0 += sy;
		}
	}
}

static void quad_bezier_sgm(int x0, int y0, int x1, int y1,
			    int xcontrol, int ycontrol,
			    void (*set_pixel)(int x, int y,
					      uint8_t i, void*),
			    void *pixel_data)
{
	int sx = x1 - xcontrol;
	int sy = y1 - ycontrol;
	long xx = x0 - xcontrol;
	long yy = y0 - ycontrol;
	double cur = xx * sy - yy * sx;        /* curvature */

	assert(xx * sx <= 0 && yy * sy <= 0);  /* Gradient sign must not change */

	if (sx * (long)sx + sy * (long)sy > xx * xx + yy * yy) {
		x1 = x0;
		x0 = sx + xcontrol;
		y1 = y0;
		y0 = sy + ycontrol;
		cur = -cur;
	}
	if (cur != 0) {
		xx += sx;
		xx *= sx = x0 < x1 ? 1 : -1;
		yy += sy;
		yy *= sy = y0 < y1 ? 1 : -1;
		double xy = 2 * xx * yy;
		xx *= xx;
		yy *= yy;
		if (cur * sx * sy < 0) {
			xx = -xx;
			yy = -yy;
			xy = -xy;
			cur = -cur;
		}
		double dx = 4.0 * sy * cur * (xcontrol - x0) + xx - xy;
		double dy = 4.0 * sx * cur * (y0 - ycontrol) + yy - xy;
		xx += xx;
		yy += yy;
		double err = dx + dy + xy;
		do {
			set_pixel(x0, y0, 255, pixel_data);
			if (x0 == x1 && y0 == y1)
				return;
			ycontrol = 2 * err < dx;
			if (2 * err > dy) {
				x0 += sx;
				dx -= xy;
				dy += yy;
				err += dy;
			}
			if (ycontrol) {
				y0 += sy;
				dy -= xy;
				dx += xx;
				err += dx;
			}
		} while (dy < 0 && dx > 0);
	}
	line(x0, y0, x1, y1, set_pixel, pixel_data);
}


static void quad_bezier(int x0, int y0, int x1, int y1,
			int xcontrol, int ycontrol, 
			void (*set_pixel)(int x, int y,
					  uint8_t i, void*),
			void *pixel_data,
			void (*qb_sgm)(int, int, int, int, int, int,
				       void (*set_pixel)
				       (int, int, uint8_t, void*),
				       void*))
{
	int x = x0 - xcontrol;
	int y = y0 - ycontrol;
	double t = x0 - 2 * xcontrol + x1;

	if ((long)x * (x1 - xcontrol) > 0) {
		if ((long)y * (y1 - ycontrol) > 0) {
			if (fabs((y0 - 2 * ycontrol + y1) / t * x) > abs(y)) {
				x0 = x1;
				x1 = x + xcontrol;
				y0 = y1;
				y1 = y + ycontrol;
			}
		}
		t = (x0 - xcontrol) / t;
		double r = (1-t) * ((1-t) * y0 + 2.0 * t * ycontrol) + POW2(t) * y1;
		t = (x0 * x1 - POW2(xcontrol)) * t / (x0 - xcontrol);
		x = floor(t + 0.5);
		y = floor(r + 0.5);
		r = (ycontrol - y0) * (t - x0) / (xcontrol - x0) + y0;
		qb_sgm(x0, y0, x, y, x, floor(r + 0.5),
		       set_pixel, pixel_data);
		r = (ycontrol - y1) * (t - x1) / (xcontrol - x1) + y1;
		x0 = x;
		xcontrol = x;
		y0 = y;
		ycontrol = floor(r + 0.5);
	}                                                 
	if ((long)(y0 - ycontrol) * (y1 - ycontrol) > 0) {
		t = y0 - 2 * ycontrol + y1;
		t = (y0 - ycontrol) / t;
		double r = (1-t) * ((1-t) * x0 + 2.0 * t * xcontrol) + POW2(t) * x1;
		t = (y0 * y1 - POW2(ycontrol)) * t / (y0 - ycontrol);
		x = floor(r + 0.5);
		y = floor(t + 0.5);
		r = (xcontrol - x0) * (t - y0) / (ycontrol - y0) + x0;
		qb_sgm(x0, y0, x, y, floor(r+0.5), y,
				set_pixel, pixel_data);
		r = (xcontrol - x1) * (t - y1) / (ycontrol - y1) + x1;
		x0 = x;
		xcontrol = floor(r + 0.5);
		y0 = y;
		ycontrol = y;
	}
	qb_sgm(x0, y0, x1, y1,
	       xcontrol, ycontrol,
	       set_pixel, pixel_data);
}

static void quad_rational_bezier_sgm
				(int x0, int y0, int x1, int y1,
				 int xcontrol, int ycontrol, float w, 
				 void (*set_pixel)(int x, int y,
						   uint8_t i, void*),
				 void *pixel_data)
{
	int sx = x1 - xcontrol;
	int sy = y1 - ycontrol;
	double dx = x0 - x1;
	double dy = y0 - y1;
	double xx = x0 - xcontrol;
	double yy = y0 - ycontrol;
	double xy = xx * sy + yy * sx;
	double cur = xx * sy - yy * sx;

	assert(xx * sx <= 0.0 && yy * sy <= 0.0);

	if (cur != 0.0 && w > 0.0) {
		if (sx * (long)sx + sy * (long)sy > xx * xx + yy * yy) {
			x1 = x0;
			x0 -= dx;
			y1 = y0;
			y0 -= dy;
			cur = -cur;
		}
		xx = 2.0 * (4.0 * w * sx * xx + dx * dx);
		yy = 2.0 * (4.0 * w * sy * yy + dy * dy);
		sx = x0 < x1 ? 1 : -1;
		sy = y0 < y1 ? 1 : -1;
		xy = -2.0 * sx * sy * (2.0 * w * xy + dx * dy);

		if (cur * sx * sy < 0.0) {
			xx = -xx;
			yy = -yy;
			xy = -xy;
			cur = -cur;
		}
		dx = 4.0 * w * (xcontrol - x0) * sy * cur + xx / 2.0 + xy;
		dy = 4.0 * w * (y0 - ycontrol) * sx * cur + yy / 2.0 + xy;

		if (w < 0.5 && (dy > xy || dx < xy)) {
			cur = (w + 1.0) / 2.0;
			w = sqrt(w);
			xy = 1.0/(w + 1.0);
			sx = floor((x0 + 2.0*w*xcontrol + x1) * xy / 2.0 + 0.5);
			sy = floor((y0 + 2.0*w*ycontrol + y1) * xy / 2.0 + 0.5);
			dx = floor((w * xcontrol + x0) * xy + 0.5);
			dy = floor((ycontrol * w + y0) * xy + 0.5);
			quad_rational_bezier_sgm(x0, y0, dx, dy, sx, sy, cur,
						 set_pixel, pixel_data);
			dx = floor((w * xcontrol + x1) * xy + 0.5);
			dy = floor((ycontrol * w + y1) * xy + 0.5);
			quad_rational_bezier_sgm(sx, sy, dx, dy, x1, y1, cur,
						 set_pixel, pixel_data);
			return;
		}
		double err = dx + dy - xy;
		do {
			set_pixel(x0, y0, 255, pixel_data);
			if (x0 == x1 && y0 == y1)
				return;
			xcontrol = 2 * err > dy;
			ycontrol = 2 * (err + yy) < -dy;
			if (2 * err < dx || ycontrol) {
				y0 += sy;
				dy += xy;
				err += dx += xx;
			}
			if (2 * err > dx || xcontrol) {
				x0 += sx;
				dx += xy;
				err += dy += yy;
			}
		} while (dy <= xy && dx >= xy);
	}
	line(x0, y0, x1, y1, set_pixel, pixel_data);
}

static void quad_rational_bezier(int x0, int y0, int x1, int y1,
				 int xcontrol, int ycontrol, float w,
				 void (*set_pixel)(int x, int y,
						   uint8_t i, void*),
				 void *pixel_data,
				 void (*qrb_sgm)(int, int, int, int,
						 int, int, float, 
						 void (*set_pixel)
						 (int, int, uint8_t, void*), 	
						 void *))
{
	int x = x0 - 2 * xcontrol + x1;
	int y = y0 - 2 * ycontrol + y1;
	double xx = x0 - xcontrol;
	double yy = y0 - ycontrol;

	if (w <= 0.0)
		w = fabs(w);

	if (xx * (x1 - xcontrol) > 0) {
		if (yy * (y1 - ycontrol) > 0) {
			if (fabs(xx * y) > fabs(yy * x)) {
				x0 = x1;
				x1 = xx + xcontrol;
				y0 = y1;
				y1 = yy + ycontrol;
			}
		}
		double t, q;
		if (x0 == x1 || w == 1.0) {
			t = (x0 - xcontrol) / (double)x;
		} else {
			q = sqrt(4.0 * POW2(w) * (x0 - xcontrol) *
				 (x1 - xcontrol) +
				 (x1 - x0) * (long)(x1 - x0));
			if (xcontrol < x0)
				q = -q;
			t = (2.0 * w * (x0 - xcontrol) - x0 + x1 + q) /
				(2.0 * (1.0 - w) * (x1 - x0));
		}
		q = 1.0 / (2.0 * t * (1.0 - t) * (w - 1.0) + 1.0);
		xx = (POW2(t) * (x0 - 2.0 * w * xcontrol + x1) +
		      2.0 * t * (w * xcontrol - x0) + x0) * q;
		yy = (POW2(t) * (y0 - 2.0 * w * ycontrol + y1) +
		      2.0 * t * (w * ycontrol - y0) + y0) * q;
		double ww = t * (w - 1.0) + 1.0;
		ww *= ww * q;
		w = ((1.0 - t) * (w - 1.0) + 1.0) * sqrt(q);
		x = floor(xx + 0.5);
		y = floor(yy + 0.5);
		yy = (xx - x0) * (ycontrol - y0) / (xcontrol - x0) + y0;
		qrb_sgm(x0, y0, x, y, 
			x, floor(yy + 0.5), ww,
			set_pixel, pixel_data);
		yy = (xx - x1) * (ycontrol - y1) / (xcontrol - x1) + y1;
		ycontrol = floor(yy + 0.5);
		x0 = xcontrol = x;
		y0 = y;
	}
	if ((y0 - ycontrol) * (long)(y1 - ycontrol) > 0) {
		double t, q;
		if (y0 == y1 || w == 1.0) {
			t = (y0 - ycontrol) / (y0 - 2.0 * ycontrol + y1);
		} else {
			q = sqrt(4.0 * POW2(w) * (y0 - ycontrol) *
				 (y1 - ycontrol) +
				 (y1 - y0) * (long)(y1 - y0));
			if (ycontrol < y0)
				q = -q;
			t = (2.0 * w * (y0 - ycontrol) - y0 + y1 + q) /
				(2.0 * (1.0 - w) * (y1 - y0));
		}
		q = 1.0 / (2.0 * t * (1.0 - t) * (w - 1.0) + 1.0);
		xx = (POW2(w) * (x0 - 2.0 * w * xcontrol + x1) +
		      2.0 * t * (w * xcontrol - x0) + x0) * q;
		yy = (POW2(w) * (y0 - 2.0 * w * ycontrol + y1) +
		      2.0 * t * (w * ycontrol - y0) + y0) * q;
		double ww = t * (w - 1.0) + 1.0;
		ww *= ww * q;
		w = ((1.0 - t) * (w - 1.0) + 1.0) * sqrt(q);
		x = floor(xx + 0.5);
		y = floor(yy + 0.5);
		xx = (xcontrol - x0) * (yy - y0) / (ycontrol - y0) + x0;
		qrb_sgm(x0, y0, x, y,
			floor(xx + 0.5), y, ww,
			set_pixel, pixel_data);
		xx = (xcontrol - x1) * (yy - y1) /
			(ycontrol - y1) + x1;
		xcontrol = floor(xx + 0.5);
		x0 = x;
		y0 = ycontrol = y;
	}
	qrb_sgm(x0, y0, x1, y1, 
		xcontrol, ycontrol, POW2(w),
		set_pixel, pixel_data);
}

static void cubic_bezier_sgm(int x0, int y0, int x1, int y1,
			     float x0_control, float y0_control,
			     float x1_control, float y1_control,
			     void (*set_pixel)(int x, int y,
					       uint8_t i, void*),
			     void *pixel_data)
{
	int leg = 1;
	int sx = x0 < x1 ? 1 : -1;
	int sy = y0 < y1 ? 1 : -1;
	float xc = -fabs(x0+x0_control-x1_control-x1);
	float xa = xc-4*sx*(x0_control-x1_control);
	float xb = sx*(x0-x0_control-x1_control+x1);
	float yc = -fabs(y0+y0_control-y1_control-y1);
	float ya = yc-4*sy*(y0_control-y1_control);
	float yb = sy*(y0-y0_control-y1_control+y1);
	double EP = 0.01;

	assert((x0_control-x0) * (x1_control-x1) < EP &&
	       ((x1-x0) * (x0_control-x1_control) < EP ||
		POW2(xb) < xa*xc + EP));

	assert((y0_control-y0) * (y1_control-y1) < EP &&
	       ((y1-y0) * (y0_control-y1_control) < EP ||
		POW2(yb) < ya*yc + EP));

	if (xa == 0 && ya == 0) {
		sx = floor((3*x0_control-x0+1)/2);
		sy = floor((3*y0_control-y0+1)/2);
		return quad_bezier_sgm(x0, y0, x1, y1, sx, sy,
				       set_pixel, pixel_data);
	}
	x0_control = (x0_control-x0) * (x0_control-x0) +
		(y0_control-y0) * (y0_control-y0) + 1;
	x1_control = (x1_control-x1) * (x1_control-x1) +
		(y1_control-y1) * (y1_control-y1) + 1;
	do {
		double ab = xa*yb - xb*ya;
		double ac = xa*yc - xc*ya;
		double bc = xb*yc - xc*yb;
		double ex = ab*(ab + ac - 3*bc) + ac*ac;
		int f = (ex > 0) ? 1 : sqrt(1 + 1024/x0_control);
		ab *= f;
		ac *= f;
		bc *= f;
		ex *= POW2(f);
		double xy = 9*(ab + ac + bc)/8;
		double cb = 8*(xa - ya);
		double dx = 27*(8*ab*(POW2(yb) - ya*yc) +
				ex*(ya + 2*yb+yc)) / 64 - POW2(ya) * (xy-ya);
		double dy = 27*(8*ab*(POW2(xb) - xa*xc) -
				ex*(xa + 2*xb+xc)) / 64 - POW2(xa) * (xy+xa);
		double xx = 3 * (3 * ab * (3 * POW2(yb) - POW2(ya) - 2*ya*yc) -
				 ya*(3*ac*(ya + yb) + ya*cb))/4;
		double yy = 3 * (3 * ab * (3 * POW2(xb) - POW2(xa) - 2*xa*xc) -
				 xa*(3*ac*(xa + xb) + xa*cb))/4;
		xy = xa * ya * (6*ab + 6*ac - 3*bc + cb);
		ac = POW2(ya);
		cb = POW2(xa);
		xy = 3*(xy + 9*f*(cb*yb*yc - xb*xc*ac) - 18*xb*yb*ab)/8;

		if (ex < 0) {
			dx = -dx;
			dy = -dy;
			xx = -xx;
			yy = -yy;
			xy = -xy;
			ac = -ac;
			cb = -cb;
		}
		ab = 6 * ya * ac;
		ac = -6 * xa * ac;
		bc = 6 * ya * cb;
		cb = -6 * xa * cb;
		dx += xy;
		ex = dx + dy;
		dy += xy;

		double *pxy = &xy;
		int fx = f;
		int fy = f;
		while (x0 != x1 && y0 != y1) {
			set_pixel(x0, y0, 255, pixel_data);
			do {
				if (dx > *pxy || dy < *pxy)
					goto exit;
				y0_control = 2*ex - dy;
				if (2*ex >= dx) {
					fx--;
					dx += xx;
					ex += dx;
					xy += ac;
					dy += xy;
					yy += bc;
					xx += ab;
				}
				if (y0_control <= 0) {
					fy--;
					dy += yy;
					ex += dy;
					xy += bc;
					dx += xy;
					xx += ac;
					yy += cb;
				}
			} while (fx > 0 && fy > 0);
			if (2*fx <= f) {
				x0 += sx;
				fx += f;
			}
			if (2*fy <= f) {
				y0 += sy;
				fy += f;
			}
			if (pxy == &xy && dx < 0 && dy > 0)
				pxy = &EP;
		}
	exit:
		xx = x0;
		x0 = x1;
		x1 = xx;
		sx = -sx;
		xb = -xb;
		
		yy = y0;
		y0 = y1;
		y1 = yy;
		sy = -sy;
		yb = -yb;
		x0_control = x1_control;
	} while (leg--);

	line(x0, y0, x1, y1, set_pixel, pixel_data);
}

static void cubic_bezier(int x0, int y0, int x1, int y1,
			 int x0_control, int y0_control,
			 int x1_control, int y1_control,
			 void (*set_pixel)(int x, int y,
					   uint8_t i, void*),
			 void *pixel_data,
			 void (*cb_sgm)(int, int, int, int,
					float, float, float, float,
					void (*set_pixel)
					(int, int, uint8_t, void*),
					void*))
{
	long xc = x0 + x0_control - x1_control - x1;
	long xa = xc - 4 * (x0_control - x1_control);
	long xb = x0 - x0_control - x1_control + x1;
	long xd = xb + 4 * (x0_control + x1_control);
	long yc = y0 + y0_control - y1_control - y1;
	long ya = yc - 4 * (y0_control - y1_control);
	long yb = y0 - y0_control - y1_control + y1;
	long yd = yb + 4 * (y0_control + y1_control);
	double t1 = xb * xb - xa * xc;
	double t[5];

	/* Sub-divide curve at gradient sign changes */
	int n = 0;
	if (xa == 0) {
		if (abs(xc) < 2*abs(xb))
			t[n++] = xc / (2.0 * xb);
	} else if (t1 > 0.0) {
		double t2 = sqrt(t1);
		t1 = (xb - t2) / xa;
		if (fabs(t1) < 1.0)
			t[n++] = t1;
		t1 = (xb+t2)/xa;
		if (fabs(t1) < 1.0)
			t[n++] = t1;
	}
	t1 = yb * yb - ya * yc;
	if (ya == 0) {
		if (abs(yc) < 2*abs(yb))
			t[n++] = yc / (2.0 * yb);
	} else if (t1 > 0.0) {
		double t2 = sqrt(t1);
		t1 = (yb-t2) / ya;
		if (fabs(t1) < 1.0)
			t[n++] = t1;
		t1 = (yb+t2) / ya;
		if (fabs(t1) < 1.0)
			t[n++] = t1;
	}
	for (int i = 1; i < n; i++) {
		if ((t1 = t[i-1]) > t[i]) {
			t[i-1] = t[i];
			t[i] = t1;
			i = 0;
		}
	}

	t1 = -1.0;
	t[n] = 1.0;
	for (int i = 0; i <= n; i++) {
		float fx0 = x0;
		float fy0 = y0;
		float fx1 = x1;
		float fy1 = y1;
		double t2 = t[i];
		float fx0_control = (t1*(t1*xb-2*xc) -
				     t2*(t1*(t1*xa-2*xb) + xc) + xd)/8 - fx0;
		float fy0_control = (t1*(t1*yb-2*yc) -
				     t2*(t1*(t1*ya-2*yb) + yc) + yd)/8 - fy0;
		float fx1_control = (t2*(t2*xb-2*xc) -
				     t1*(t2*(t2*xa-2*xb) + xc) + xd)/8 - fx0;
		float fy1_control = (t2*(t2*yb-2*yc) -
				     t1*(t2*(t2*ya-2*yb) + yc) + yd)/8 - fy0;

		fx0 -= fx1 = (t2*(t2*(3*xb - t2*xa) - 3*xc) + xd) / 8;
		fy0 -= fy1 = (t2*(t2*(3*yb - t2*ya) - 3*yc) + yd) / 8;
		x1 = floor(fx1 + 0.5);
		y1 = floor(fy1 + 0.5);
		if (fx0 != 0.0) {
			fx0 = (x0 - x1) / fx0;
			fx0_control *= fx0;
			fx1_control *= fx0;
		}
		if (fy0 != 0.0) {
			fy0 = (y0 - y1) / fy0;
			fy0_control *= fy0;
			fy1_control *= fy0;
		}
		if (x0 != x1 || y0 != y1)
			cb_sgm(x0, y0, x1, y1,
			       x0 + fx0_control,
			       y0 + fy0_control,
			       x0 + fx1_control,
			       y0 + fy1_control,
			       set_pixel,
			       pixel_data);
		x0 = x1;
		y0 = y1;
		fx0 = fx1;
		fy0 = fy1;
		t1 = t2;
	}
}

static void line_aa(int x0, int y0, int x1, int y1,
		    void (*set_pixel)(int x, int y,
				      uint8_t i, void*),
		    void *pixel_data)
{
	int sx = (x0 < x1)? 1 : -1;
	int sy = (y0 < y1)? 1 : -1;
	long dx = abs(x1 - x0);
	long dy = abs(y1 - y0);
	long err = POW2(dx) + POW2(dy);
	long e2 = (err == 0)? 1 : 0xffff7fl / sqrt(err);

	dx *= e2;
	dy *= e2;
	err = dx - dy;
	while (true) {
		set_pixel(x0, y0, 255 - (abs(err-dx+dy)>>16), pixel_data);
		e2 = err;
		int x2 = x0;
		if (2*e2 >= -dx) {
			if (x0 == x1)
				break;
			if (e2+dy < 0xff0000l)
				set_pixel(x0, y0 + sy,
					  255 - ((e2+dy)>>16),
					  pixel_data);
			err -= dy;
			x0 += sx; 
		} 
		if (2*e2 <= dy) {
			if (y0 == y1)
				break;
			if (dx-e2 < 0xff0000l)
				set_pixel(x2+sx, y0, 255 -
					  ((dx-e2)>>16),
					  pixel_data);
			err += dx;
			y0 += sy; 
		}
	}
}

static void quad_bezier_sgm_aa(int x0, int y0, int x1, int y1,
			       int xcontrol, int ycontrol,
			       void (*set_pixel)(int x, int y,
						 uint8_t i, void*),
			       void *pixel_data)
{
	int sx = x1 - xcontrol;
	int sy = y1 - ycontrol;
	long xx = x0 - xcontrol;
	long yy = y0 - ycontrol;
	double cur = xx*sy - yy*sx;

	assert(xx*sx <= 0 && yy*sy <= 0);

	if (sx*(long)sx + sy*(long)sy > POW2(xx) + POW2(yy)) {
		x1 = x0;
		x0 = sx + xcontrol;
		y1 = y0;
		y0 = sy + ycontrol;
		cur = -cur;
	}
	if (cur != 0) {
		xx += sx;
		sx = (x0 < x1)? 1 : -1;
		xx *= sx;
		yy += sy;
		sy = (y0 < y1)? 1 : -1;
		yy *= sy;
		long xy = 2*xx*yy;
		xx *= xx;
		yy *= yy;
		if (cur*sx*sy < 0) {
			xx = -xx;
			yy = -yy;
			xy = -xy;
			cur = -cur;
		}
		double dx = 4.0*sy*(xcontrol-x0)*cur + xx - xy;
		double dy = 4.0*sx*(y0-ycontrol)*cur + yy - xy;
		xx += xx;
		yy += yy;
		double err = dx + dy + xy;
		do {
			cur = fmin(dx + xy, -xy - dy);
			double ed = fmax(dx + xy, -xy - dy);
			ed += 2*ed* POW2(cur) / (POW2(2*ed) + POW2(cur));
			set_pixel(x0, y0, 255 * (1.0 - fabs(err-dx-dy-xy)/ed),
				  pixel_data);
			if (x0 == x1 || y0 == y1)
				break;
			xcontrol = x0;
			cur = dx - err;
			ycontrol = 2*err+dy < 0;
			if (2*err+dx > 0) {
				if (err-dy < ed)
					set_pixel(x0, y0 + sy,
						  255 * (1.0 - fabs(err-dy)/ed),
						  pixel_data);
				x0 += sx;
				dx -= xy;
				dy += yy;
				err += dy;
			}
			if (ycontrol) {
				if (cur < ed)
					set_pixel(xcontrol+sx, y0,
						  255 * (1.0 - fabs(cur)/ed),
						  pixel_data);
				y0 += sy;
				dy -= xy;
				dx += xx;
				err += dx;
			}
		} while (dy < dx);
	}
	line_aa(x0, y0, x1, y1, set_pixel, pixel_data);
}

static void quad_rational_bezier_sgm_aa
				(int x0, int y0, int x1, int y1,
				 int xcontrol, int ycontrol, float w, 
				 void (*set_pixel)(int x, int y,
						   uint8_t i, void*),
				 void *pixel_data)
{
	int sx = x1 - xcontrol;
	int sy = y1 - ycontrol;
	double dx = x0 - x1;
	double dy = y0 - y1;
	double xx = x0 - xcontrol;
	double yy = y0 - ycontrol;
	double xy = xx*sy + yy*sx;
	double cur = xx*sy - yy*sx;

	assert(xx*sx <= 0.0 && yy*sy <= 0.0);

	if (cur != 0.0 && w > 0.0) {
		if (sx*(long)sx + sy*(long)sy > POW2(xx) + POW2(yy)) {
			x1 = x0;
			x0 -= dx;
			y1 = y0;
			y0 -= dy;
			cur = -cur;
		}
		xx = 2.0 * (4.0*w*sx*xx + POW2(dx));
		yy = 2.0 * (4.0*w*sy*yy + POW2(dy));
		sx = (x0 < x1)? 1 : -1;
		sy = (y0 < y1)? 1 : -1;
		xy = -2.0*sx*sy * (2.0*w*xy + dx*dy);

		if (cur*sx*sy < 0) {
			xx = -xx;
			yy = -yy;
			cur = -cur;
			xy = -xy;
		}
		dx = 4.0*w * (xcontrol - x0) * sy*cur + xx/2.0 + xy;
		dy = 4.0*w * (y0 - ycontrol) * sx*cur + yy/2.0 + xy;

		if (w < 0.5 && dy > dx) {
			cur = (w+1.0)/2.0;
			w = sqrt(w);
			xy = 1.0/(w+1.0);
			sx = floor((x0 + 2.0*w*xcontrol + x1) * xy/2.0 + 0.5);
			sy = floor((y0 + 2.0*w*ycontrol + y1) * xy/2.0 + 0.5);
			dx = floor((w*xcontrol + x0) * xy + 0.5);
			dy = floor((w*ycontrol + y0) * xy + 0.5);
			quad_rational_bezier_sgm_aa(x0, y0, sx,sy, dx,dy, cur,
						    set_pixel, pixel_data);
			dx = floor((w*xcontrol + x1) * xy + 0.5);
			dy = floor((w*ycontrol + y1) * xy + 0.5);
			quad_rational_bezier_sgm_aa(sx,sy, x1,y1,
						    dx, dy, cur,
						    set_pixel,
						    pixel_data);
			return;
		}
		double err = dx + dy - xy;
		do {
			cur = fmin(dx-xy, xy-dy);
			double ed = fmax(dx-xy, xy-dy);
			ed += 2*ed*POW2(cur)/(POW2(2.0*ed) + POW2(cur));
			xcontrol = 255 * fabs(err - dx - dy + xy)/ed;
			if (xcontrol < 256)
				set_pixel(x0, y0, 255 - xcontrol, pixel_data);
			bool f = (2*err + dy < 0);
			if (f) {
				if (y0 == y1)
					return;
				if (dx-err < ed)
					set_pixel(x0 + sx, y0,
						  255 * (1.0 - fabs(dx-err)/ed),
						  pixel_data);
			}
			if (2*err+dx > 0) {
				if (x0 == x1)
					return;
				if (err-dy < ed)
					set_pixel(x0, y0 + sy,
						  255 * (1.0 - fabs(err-dy)/ed),
						  pixel_data);
				x0 += sx;
				dx += xy;
				dy += yy;
				err += dy;
			}
			if (f) {
				y0 += sy;
				dy += xy;
				dx += xx;
				err += dx;
			}
		} while (dy < dx);
	}
	line_aa(x0, y0, x1, y1, set_pixel, pixel_data);
}

static void cubic_bezier_sgm_aa(int x0, int y0, int x1, int y1,
				float x0_control, float y0_control,
				float x1_control, float y1_control,
				void (*set_pixel)(int x, int y,
						  uint8_t i, void*),
				void *pixel_data)
{
	int leg = 1;
	int sx = (x0 < x1)? 1 : -1;
	int sy = (y0 < y1)? 1 : -1;
	float xc = -fabs(x0 + x0_control - x1_control - x1);
	float xa = xc - 4*sx*(x0_control - x1_control);
	float xb = sx * (x0 - x0_control - x1_control + x1);
	float yc = -fabs(y0 + y0_control - y1_control - y1);
	float ya = yc - 4*sy*(y0_control - y1_control);
	float yb = sy*(y0 - y0_control - y1_control + y1);
	double EP = 0.01;

	assert((x0_control-x0)*(x1_control-x1) < EP &&
	       ((x1-x0)*(x0_control-x1_control) < EP || xb*xb < xa*xc+EP));
	assert((y0_control-y0)*(y1_control-y1) < EP &&
	       ((y1-y0)*(y0_control-y1_control) < EP || yb*yb < ya*yc+EP));

	if (xa == 0 && ya == 0) {
		sx = floor((3*x0_control - x0 + 1) / 2);
		sy = floor((3*y0_control - y0 + 1) / 2);
		quad_bezier_sgm_aa(x0,y0, x1,y1, sx,sy,
				   set_pixel, pixel_data);
		return;
	}
	x0_control = (x0_control-x0)*(x0_control-x0) + 
		(y0_control-y0)*(y0_control-y0) + 1;
	x1_control = (x1_control-x1)*(x1_control-x1) +
		(y1_control-y1)*(y1_control-y1) + 1;
	do {
		double ab = xa*yb - xb*ya;
		double ac = xa*yc - xc*ya;
		double bc = xb*yc - xc*yb;
		double ip = 4*ab*bc - POW2(ac);
		double ex = ab*(ab+ac-3*bc) + POW2(ac);
		int f = (ex > 0)? 1 : sqrt(1+1024/x0_control);
		ab *= f;
		ac *= f;
		bc *= f;
		ex *= POW2(f);
		double xy = 9*(ab + ac + bc)/8;
		double ba = 8*(xa - ya);
		double dx = 27*(8*ab*(yb*yb-ya*yc) +
				ex*(ya+2*yb+yc))/64-ya*ya*(xy-ya);
		double dy = 27*(8*ab*(xb*xb-xa*xc) -
				ex*(xa+2*xb+xc))/64-xa*xa*(xy+xa);

		double xx = 3*(3*ab*(3*yb*yb-ya*ya-2*ya*yc) -
			       ya*(3*ac*(ya+yb)+ya*ba))/4;
		double yy = 3*(3*ab*(3*xb*xb-xa*xa-2*xa*xc) -
			       xa*(3*ac*(xa+xb)+xa*ba))/4;
		xy = xa*ya*(6*ab+6*ac-3*bc+ba);
		ac = POW2(ya);
		ba = POW2(xa);
		xy = 3*(xy+9*f*(ba*yb*yc-xb*xc*ac)-18*xb*yb*ab)/8;

		if (ex < 0) {
			dx = -dx;
			dy = -dy;
			xx = -xx;
			yy = -yy;
			xy = -xy;
			ac = -ac;
			ba = -ba;
		}
		ab = 6*ya*ac;
		ac = -6*xa*ac;
		bc = 6*ya*ba;
		ba = -6*xa*ba;
		dx += xy;
		ex = dx + dy;
		dy += xy;

		int fx = f;
		int fy = f;
		double ed, px, py;
		while (x0 != x1 && y0 != y1) {
			y0_control = fmin(fabs(xy-dx), fabs(dy-xy));
			ed = fmax(fabs(xy-dx), fabs(dy-xy));
			ed = f*(ed + 2*ed*y0_control*y0_control /
				(4*ed*ed+y0_control*y0_control));
			y0_control = 255 * fabs(ex-(f-fx+1)*dx -
						(f-fy+1)*dy+f*xy)/ed;
			if (y0_control < 256)
				set_pixel(x0, y0, 255 - y0_control, pixel_data);
			px = fabs(ex-(f-fx+1)*dx + (fy-1)*dy);
			py = fabs(ex+(fx-1)*dx - (f-fy+1)*dy);
			y1_control = y0;
			do {
				if (ip >= -EP)
					if (dx + xx > xy || dy + yy < xy)
						goto exit;
				y0_control = 2*ex + dx;
				if (2*ex+dy > 0) {
					fx -= 1;
					dx += xx;
					ex += dx;
					xy += ac;
					dy += xy;
					yy += bc;
					xx += ab;
				} else if (y0_control > 0) {
					goto exit;
				}
				if (y0_control <= 0) {
					fy--;
					dy += yy;
					ex += dy;
					xy += bc;
					dx += xy;
					xx += ac;
					yy += ba;
				}
			} while (fx > 0 && fy > 0);
			if (2*fy <= f) {
				if (py < ed)
					set_pixel(x0 + sx, y0,
						  255 * (1.0 - py/ed), 
						  pixel_data);
				y0 += sy;
				fy += f;
			}
			if (2*fx <= f) {
				if (px < ed)
					set_pixel(x0, y1_control + sy,
						  255 * (1.0 - px/ed),
						  pixel_data);
				x0 += sx;
				fx += f;
			}
		}
		break;
	exit:
		if (2*ex < dy && 2*fy <= f+2) {
			if (py < ed)
				set_pixel(x0 + sx, y0,
					  255 * (1.0 - py/ed),
					  pixel_data);
			y0 += sy;
		}
		if (2*ex > dx && 2*fx <= f+2) {
			if (px < ed)
				set_pixel(x0, y1_control + sy,
					  255 * (1.0 - px/ed),
					  pixel_data);
			x0 += sx;
		}
		xx = x0;
		x0 = x1;
		x1 = xx;
		sx = -sx;
		xb = -xb;
		yy = y0;
		y0 = y1;
		y1 = yy;
		sy = -sy;
		yb = -yb;
		x0_control = x1_control;
	} while (leg--);

	line_aa(x0, y0, x1, y1, set_pixel, pixel_data);
}

void nb_graphics_rasterizer_fill(int x, int y, uint8_t i,
				 void (*set_pixel)(int x, int y,
						   uint8_t i, void*),
				 bool (*pixel_is_empty)(int x, int y,
							const void*),
				 uint32_t width, uint32_t height,
				 void *pixel_data)
/* 4-way scanline */
{
	if (!pixel_is_empty(x, y, pixel_data))
		goto EXIT;

	pix_stack_t *pix_stack = alloca(sizeof(*pix_stack));
	init_stack(pix_stack, width);

	push(pix_stack, x, y);

	while (pop(pix_stack, &x, &y)) {
		int x1 = x;
		while (x1 >= 0 && pixel_is_empty(x1, y, pixel_data))
			x1--;
		x1++;

		bool span_above = false;
		bool span_below = false;
		while (x1 < width && pixel_is_empty(x1, y, pixel_data)) {
			set_pixel(x1, y, i, pixel_data);
			if (y > 0) {
				bool pix_is_empty =
					pixel_is_empty(x1, y - 1, pixel_data);
				if (!span_above && pix_is_empty) {
					push(pix_stack, x1, y - 1);
					span_above = true;
				} else if (span_above && !pix_is_empty) {
					span_above = false;
				}
			}
			if (y < height - 1) {
				bool pix_is_empty =
					pixel_is_empty(x1, y + 1, pixel_data);
				if (!span_below && pix_is_empty) {
					push(pix_stack, x1, y + 1);
					span_below = true;
				} else if (span_below && !pix_is_empty) {
					span_below = false;
				}
			}
			x1++;
		}
	}
	finish_stack(pix_stack);
 EXIT:
	return;
}

static void init_stack(pix_stack_t *ps, int32_t img_width)
{
	ps->w = img_width;
	ps->N = 0;
	ps->N_alloc = 0;
	memset(ps->static_mem, 0,
	       PIX_STACK_STATIC_SIZE * sizeof(*(ps->static_mem)));
	ps->dynamic_mem = NULL;
}

static bool pop(pix_stack_t *ps, int *x, int *y)
{
	bool out = false;
	if (0 < ps->N) {
		uint32_t p;
		if (ps->N <= PIX_STACK_STATIC_SIZE) {
			uint32_t last = ps->N - 1;
			p = ps->static_mem[last];
		} else {
			uint32_t last = ps->N - PIX_STACK_STATIC_SIZE - 1;
			p = ps->dynamic_mem[last];
		}
		ps->N -= 1;
		*x = p % ps->w;
		*y = p / ps->w;
		out = true;
	}
	return out;
}

static void push(pix_stack_t *ps, int x, int y)
{
	uint32_t p = y * ps->w + x;
	if (ps->N < PIX_STACK_STATIC_SIZE) {
		ps->static_mem[ps->N] = p;
	} else {
		uint32_t N = ps->N - PIX_STACK_STATIC_SIZE;
		if (ps->N_alloc <= N) {
			ps->N_alloc += PIX_STACK_DYNAMIC_INCREMENTS;
			uint32_t memsize = ps->N_alloc *
				sizeof(*(ps->dynamic_mem));
			ps->dynamic_mem =
				realloc(ps->dynamic_mem, memsize);
		}
		ps->dynamic_mem[N] = p;
	}
	ps->N += 1;
}

static void finish_stack(pix_stack_t *ps)
{
	if (ps->N_alloc > 0)
		free(ps->dynamic_mem);
}
