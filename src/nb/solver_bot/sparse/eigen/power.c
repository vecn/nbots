#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "nb/memory_bot.h"
#include "nb/solver_bot/sparse/sparse.h"
#include "nb/solver_bot/sparse/solvers/gauss_seidel.h"

#include "../sparse_struct.h"

#define POW2(a) ((a)*(a))
