#!/bin/bash

export LD_LIBRARY_PATH=core/build/libs/nbots/shared/release
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:core/build/libs/nbots_cairo/shared/release #TEMPORAL
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:core/build/libs/nb_test_add/shared/release

echo "Unit Tests"
LAUNCHER=core/build/exe/launch_unit_test/release/launch_unit_test
MODULES_DIR=core/build/libs
FLAVOR_DIR=shared/release

echo "  Container Bot"
$LAUNCHER $MODULES_DIR/container_array_UT/$FLAVOR_DIR/libcontainer_array_UT.so
$LAUNCHER $MODULES_DIR/container_QUEUE_UT/$FLAVOR_DIR/libcontainer_QUEUE_UT.so
$LAUNCHER $MODULES_DIR/container_STACK_UT/$FLAVOR_DIR/libcontainer_STACK_UT.so
$LAUNCHER $MODULES_DIR/container_SORTED_UT/$FLAVOR_DIR/libcontainer_SORTED_UT.so
$LAUNCHER $MODULES_DIR/container_HEAP_UT/$FLAVOR_DIR/libcontainer_HEAP_UT.so
$LAUNCHER $MODULES_DIR/container_HASH_UT/$FLAVOR_DIR/libcontainer_HASH_UT.so
$LAUNCHER $MODULES_DIR/iterator_QUEUE_UT/$FLAVOR_DIR/libiterator_QUEUE_UT.so
$LAUNCHER $MODULES_DIR/iterator_STACK_UT/$FLAVOR_DIR/libiterator_STACK_UT.so
$LAUNCHER $MODULES_DIR/iterator_SORTED_UT/$FLAVOR_DIR/libiterator_SORTED_UT.so
$LAUNCHER $MODULES_DIR/iterator_HEAP_UT/$FLAVOR_DIR/libiterator_HEAP_UT.so
$LAUNCHER $MODULES_DIR/iterator_HASH_UT/$FLAVOR_DIR/libiterator_HASH_UT.so

echo "  Geometric Bot"
$LAUNCHER $MODULES_DIR/point2D_UT/$FLAVOR_DIR/libpoint2D_UT.so
$LAUNCHER $MODULES_DIR/utils2D_UT/$FLAVOR_DIR/libutils2D_UT.so
$LAUNCHER $MODULES_DIR/knn_bins2D_UT/$FLAVOR_DIR/libknn_bins2D_UT.so
$LAUNCHER $MODULES_DIR/knn_bins2D_iterator_UT/$FLAVOR_DIR/libknn_bins2D_iterator_UT.so
$LAUNCHER $MODULES_DIR/mesh_dewall_UT/$FLAVOR_DIR/libmesh_dewall_UT.so
$LAUNCHER $MODULES_DIR/mesh_constrained_delaunay_UT/$FLAVOR_DIR/libmesh_constrained_delaunay_UT.so
$LAUNCHER $MODULES_DIR/model_verifier_UT/$FLAVOR_DIR/libmodel_verifier_UT.so
$LAUNCHER $MODULES_DIR/model_clipper_UT/$FLAVOR_DIR/libmodel_clipper_UT.so
$LAUNCHER $MODULES_DIR/mesh_mesh2D_UT/$FLAVOR_DIR/libmesh_mesh2D_UT.so
$LAUNCHER $MODULES_DIR/mesh_image_density_UT/$FLAVOR_DIR/libmesh_image_density_UT.so

echo "  PDE Bot"
$LAUNCHER $MODULES_DIR/pde_fem_elasticity2D_UT/$FLAVOR_DIR/libpde_fem_elasticity2D_UT.so
