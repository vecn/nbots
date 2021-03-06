#ifndef __NB_GEOMETRIC_BOT_H__
#define __NB_GEOMETRIC_BOT_H__

#include "nb/geometric_bot/point2D.h"
#include "nb/geometric_bot/utils2D.h"
#include "nb/geometric_bot/knn/bins2D.h"
#include "nb/geometric_bot/knn/bins2D_iterator.h"
#include "nb/geometric_bot/model/model2D_struct.h"
#include "nb/geometric_bot/model/model2D.h"
#include "nb/geometric_bot/model/modules2D/drawing.h"
#include "nb/geometric_bot/model/modules2D/verifier.h"
#include "nb/geometric_bot/model/modules2D/regularizer.h"
#include "nb/geometric_bot/model/modules2D/simplifier.h"
#include "nb/geometric_bot/model/modules2D/clipper.h"
#include "nb/geometric_bot/model/model3D.h"
#include "nb/geometric_bot/mesh/dewall.h"
#include "nb/geometric_bot/mesh/alpha_shape.h"
#include "nb/geometric_bot/mesh/constrained_delaunay.h"
#include "nb/geometric_bot/mesh/tessellator2D.h"
#include "nb/geometric_bot/mesh/modules2D/image_density.h"
#include "nb/geometric_bot/mesh/modules2D/graph_generator.h"
#include "nb/geometric_bot/mesh/modules2D/area_analizer.h"
#include "nb/geometric_bot/mesh/modules2D/extra_fem.h"
#include "nb/geometric_bot/mesh/modules2D/drawing.h"
#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/geometric_bot/mesh/mesh2D_file_format_vtk.h"
#include "nb/geometric_bot/mesh/mesh2D/file_format_nbt.h"
#include "nb/geometric_bot/mesh/mesh2D/info.h"
#include "nb/geometric_bot/mesh/mesh2D/draw.h"
#include "nb/geometric_bot/mesh/mesh2D/field_io.h"
#include "nb/geometric_bot/mesh/mesh2D/elements2D/trg_exporter.h"
#include "nb/geometric_bot/mesh/mesh2D/elements2D/msh3trg.h"
#include "nb/geometric_bot/mesh/mesh2D/elements2D/msh3trg_draw.h"
#include "nb/geometric_bot/mesh/mesh2D/elements2D/mshquad.h"
#include "nb/geometric_bot/mesh/mesh2D/elements2D/mshquad_draw.h"
#include "nb/geometric_bot/mesh/mesh2D/elements2D/mshpoly.h"
#include "nb/geometric_bot/mesh/mesh2D/elements2D/mshpoly_draw.h"
#include "nb/geometric_bot/mesh/mesh2D/elements2D/mshpack.h"
#include "nb/geometric_bot/mesh/mesh2D/elements2D/mshpack_draw.h"
#include "nb/geometric_bot/mesh/tessellator3D.h"

#endif
