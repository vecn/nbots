#!/bin/bash

make

echo "Unit Tests"
LAUNCHER=tests/launch_unit_tests

echo "  Container Bot"
CONTAINER_BOT_TESTS_DIR=tests/vcn/container_bot
$LAUNCHER $CONTAINER_BOT_TESTS_DIR/libarray_UT.so
$LAUNCHER $CONTAINER_BOT_TESTS_DIR/libcontainer_QUEUE_UT.so
$LAUNCHER $CONTAINER_BOT_TESTS_DIR/libcontainer_STACK_UT.so
$LAUNCHER $CONTAINER_BOT_TESTS_DIR/libcontainer_SORTED_UT.so
$LAUNCHER $CONTAINER_BOT_TESTS_DIR/libcontainer_HEAP_UT.so
$LAUNCHER $CONTAINER_BOT_TESTS_DIR/libcontainer_HASH_UT.so
$LAUNCHER $CONTAINER_BOT_TESTS_DIR/libiterator_QUEUE_UT.so
$LAUNCHER $CONTAINER_BOT_TESTS_DIR/libiterator_STACK_UT.so
$LAUNCHER $CONTAINER_BOT_TESTS_DIR/libiterator_SORTED_UT.so
$LAUNCHER $CONTAINER_BOT_TESTS_DIR/libiterator_HEAP_UT.so
$LAUNCHER $CONTAINER_BOT_TESTS_DIR/libiterator_HASH_UT.so

echo "  Geometric Bot"
GEOMETRIC_BOT_TESTS_DIR=tests/vcn/geometric_bot
$LAUNCHER $GEOMETRIC_BOT_TESTS_DIR/libpoint2D_UT.so
$LAUNCHER $GEOMETRIC_BOT_TESTS_DIR/libutils2D_UT.so
$LAUNCHER $GEOMETRIC_BOT_TESTS_DIR/libbins2D_UT.so
$LAUNCHER $GEOMETRIC_BOT_TESTS_DIR/libbins2D_iterator_UT.so
$LAUNCHER $GEOMETRIC_BOT_TESTS_DIR/libdewall_UT.so
$LAUNCHER $GEOMETRIC_BOT_TESTS_DIR/libconstrained_delaunay_UT.so
$LAUNCHER $GEOMETRIC_BOT_TESTS_DIR/libmesh2D_UT.so

echo "  Java binding for Geometric Bot"
JAVA_TEST_DIR=tests/java
JAVA_CP=src/java/VCNGeometricBot.jar:$JAVA_TEST_DIR/TestGeometricBot.jar
JAVA_LP=src:src/jni
java -Djava.library.path=$JAVA_LP -cp $JAVA_CP TestGeometricBot

