#!/bin/bash

export LD_LIBRARY_PATH=build/libs/nbots_all/shared/release
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
GEOMETRIC_BOT_TESTS_DIR=libs/*/shared/release
$LAUNCHER $MODULES_DIR/point2D_UT/$FLAVOR_DIR/libpoint2D_UT.so
#$LAUNCHER $MODULES_DIR/utils2D/$FLAVOR_DIR/libutils2D_UT.so
#$LAUNCHER $MODULES_DIR/bins2D/$FLAVOR_DIR/libbins2D_UT.so
#$LAUNCHER $MODULES_DIR/bins2D_iterator/$FLAVOR_DIR/libbins2D_iterator_UT.so
#$LAUNCHER $MODULES_DIR//$FLAVOR_DIR/libdewall_UT.so
#$LAUNCHER $MODULES_DIR//$FLAVOR_DIR/libconstrained_delaunay_UT.so
#$LAUNCHER $MODULES_DIR//$FLAVOR_DIR/libverifier_UT.so
#$LAUNCHER $MODULES_DIR//$FLAVOR_DIR/libclipper_UT.so
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
