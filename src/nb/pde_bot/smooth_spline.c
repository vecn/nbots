#include <math.h>

#include "nb/pde_bot.h"

#define POW2(a) ((a)*(a))
#define POW3(a) ((a)*(a)*(a))

double nb_smooth_spline(double x, uint8_t smooth)
{
	double spline;
	switch (smooth) {
	case 0:
		spline = x;
		break;
	case 1:
		spline = 3 * POW2(x) - 2 * POW3(x);
		break;
	case 2:
		spline = 10 * POW3(x) - 15 * pow(x, 4) + 6 * pow(x, 5);
		break;
	case 3:
		spline = 35 * pow(x, 4) - 84 * pow(x, 5) +
			70 * pow(x, 6) - 20 * pow(x, 7);
		break;
	case 4:
		spline = 126 * pow(x, 5) - 420 * pow(x, 6) +
			540 * pow(x, 7) - 315 * pow(x, 8) + 70 * pow(x, 9);
		break;
	case 5:
		spline = 462 * pow(x, 6) - 1980 * pow(x, 7) + 3465 * pow(x, 8) -
			3080 * pow(x, 9) + 1386 * pow(x, 10) - 252 * pow(x, 11);
		break;
	default:
		spline = 1716 * pow(x, 7) - 9009 * pow(x, 8) +
			20020 * pow(x, 9) - 24024 * pow(x, 10) +
			16380 * pow(x, 11) - 6006 * pow(x, 12) +
			924 * pow(x, 13);
	}
	return spline;
}
