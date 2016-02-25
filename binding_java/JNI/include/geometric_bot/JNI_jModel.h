#ifndef __JNI_NB_GEOMETRIC_BOT_JNIMODEL_H__
#define __JNI_NB_GEOMETRIC_BOT_JNIMODEL_H__

#include <jni.h>

JNIEXPORT jobject JNICALL
Java_nb_geometricBot_Model_combine(JNIEnv *env, jclass class,
				   jobject jModelA, jobject jModelB);

JNIEXPORT jobject JNICALL
Java_nb_geometricBot_Model_unify(JNIEnv *env, jclass class,
				 jobject jModelA, jobject jModelB);

JNIEXPORT jobject JNICALL
Java_nb_geometricBot_Model_intersect(JNIEnv *env, jclass class,
				     jobject jModelA, jobject jModelB);

JNIEXPORT jobject JNICALL
Java_nb_geometricBot_Model_substract(JNIEnv *env, jclass class,
				     jobject jModelA, jobject jModelB);

JNIEXPORT jobject JNICALL
Java_nb_geometricBot_Model_difference(JNIEnv *env, jclass class,
				      jobject jModelA, jobject jModelB);

JNIEXPORT jobject JNICALL
Java_nb_geometricBot_Model_verify(JNIEnv *env, jobject jModel);

JNIEXPORT jboolean JNICALL
Java_nb_geometricBot_Model_isPointInside(JNIEnv *env, jobject jModel,
					 jdouble x, jdouble y);

#endif
