#!/bin/bash

export LD_LIBRARY_PATH=build/libs/nbots_all/shared/release
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:build/libs/nbots_all_cairo/shared/release #TEMPORAL
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:build/libs/nb_test_add/shared/release

echo "Unit Tests"
LAUNCHER=build/exe/launch_unit_test/release/launch_unit_test
MODULES_DIR=build/libs
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
#$LAUNCHER $MODULES_DIR//$FLAVOR_DIR/libmesh2D_UT.so
#$LAUNCHER $MODULES_DIR//$FLAVOR_DIR/libimage_density_UT.so

#echo "  PDE Bot"
#GEOMETRIC_BOT_TESTS_DIR=libs/*/shared/release
#$LAUNCHER $GEOMETRIC_BOT_TESTS_DIR/finite_element/libstatic_elasticity2D_UT.so
#
#echo "  Java binding for Geometric Bot"
#JAVA_DIR=binding_java
#JAVA_TEST_DIR=$JAVA_DIR/tests
#JAVA_CP=$JAVA_DIR/java/GeometricBot.jar:$JAVA_DIR/java/PdeBot.jar
#JAVA_LP=$JAVA_DIR/JNI
#TEST_JARS=$JAVA_CP:$JAVA_TEST_DIR/TestModel.jar
#java -Djava.library.path=$JAVA_LP -cp $TEST_JARS TestModel
#TEST_JARS=$JAVA_CP:$JAVA_TEST_DIR/TestGeometricBot.jar
#java -Djava.library.path=$JAVA_LP -cp $TEST_JARS TestGeometricBot
#TEST_JARS=$JAVA_CP:$JAVA_TEST_DIR/TestPdeBot.jar
#java -Djava.library.path=$JAVA_LP -cp $TEST_JARS TestPdeBot
