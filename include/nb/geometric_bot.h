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
#include "nb/geometric_bot/mesh/partition.h"
#include "nb/geometric_bot/mesh/partition/info.h"
#include "nb/geometric_bot/mesh/partition/draw.h"
#include "nb/geometric_bot/mesh/partition/elements2D/trg_exporter.h"
#include "nb/geometric_bot/mesh/partition/elements2D/msh3trg.h"
#include "nb/geometric_bot/mesh/partition/elements2D/msh3trg_draw.h"
#include "nb/geometric_bot/mesh/partition/elements2D/mshquad.h"
#include "nb/geometric_bot/mesh/partition/elements2D/mshquad_draw.h"
#include "nb/geometric_bot/mesh/partition/elements2D/mshpoly.h"
#include "nb/geometric_bot/mesh/partition/elements2D/mshpoly_draw.h"
#include "nb/geometric_bot/mesh/partition/elements2D/mshpack.h"
#include "nb/geometric_bot/mesh/partition/elements2D/mshpack_draw.h"

#endif
