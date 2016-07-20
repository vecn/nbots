/**
 * @file geometric_bot.h
 * @brief The geometric bot is a program used to generate tesselations, such
 * as meshes suited for numerical interpolation in a wide range of engineering
 * applications, such as computer graphics, finite element, computer vision
 * and machine learning.
 * @author Victor Eduardo Cardoso Nungaray
 * @n victorc@@cimat.mx
 * @n <a href="https://twitter.com/victore_cardoso"> @@victore_cardoso </a>
 * @date 10 August 2015
 *
 * @mainpage Geometric Bot
 * A geometric tool for numerical analysis. Sample image:
 * @image html mesh_eye.png
 */

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
#include "nb/geometric_bot/mesh/dewall.h"
#include "nb/geometric_bot/mesh/constrained_delaunay.h"
#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/geometric_bot/mesh/modules2D/image_density.h"
#include "nb/geometric_bot/mesh/modules2D/graph_generator.h"
#include "nb/geometric_bot/mesh/modules2D/area_analizer.h"
#include "nb/geometric_bot/mesh/modules2D/extra_fem.h"
#include "nb/geometric_bot/mesh/modules2D/drawing.h"
#include "nb/geometric_bot/mesh/elements2D/trg_exporter.h"
#include "nb/geometric_bot/mesh/elements2D/triangles.h"
#include "nb/geometric_bot/mesh/elements2D/triangles_struct.h"
#include "nb/geometric_bot/mesh/elements2D/quad.h"
#include "nb/geometric_bot/mesh/elements2D/quad_struct.h"
#include "nb/geometric_bot/mesh/elements2D/polygons.h"
#include "nb/geometric_bot/mesh/elements2D/polygons_struct.h"
#include "nb/geometric_bot/mesh/elements2D/disks.h"
#include "nb/geometric_bot/mesh/elements2D/disks_struct.h"

#endif
