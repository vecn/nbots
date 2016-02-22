#ifndef __JNI_NB_GEOMETRIC_BOT_JNIMODEL_H__
#define __JNI_NB_GEOMETRIC_BOT_JNIMODEL_H__

#include <jni.h>

/*
 * Class:     nb_geometricBot_Model
 * Method:    combine
 * Signature: (Lnb/geometricBot/Model;Lnb/geometricBot/Model;)Lnb/geometricBot/Model;
 */
JNIEXPORT jobject JNICALL
Java_nb_geometricBot_Model_combine(JNIEnv *env, jclass class,
				   jobject jModelA, jobject jModelB);

/*
 * Class:     nb_geometricBot_Model
 * Method:    unify
 * Signature: (Lnb/geometricBot/Model;Lnb/geometricBot/Model;)Lnb/geometricBot/Model;
 */
JNIEXPORT jobject JNICALL
Java_nb_geometricBot_Model_unify(JNIEnv *env, jclass class,
				 jobject jModelA, jobject jModelB);

/*
 * Class:     nb_geometricBot_Model
 * Method:    intersect
 * Signature: (Lnb/geometricBot/Model;Lnb/geometricBot/Model;)Lnb/geometricBot/Model;
 */
JNIEXPORT jobject JNICALL
Java_nb_geometricBot_Model_intersect(JNIEnv *env, jclass class,
				     jobject jModelA, jobject jModelB);

/*
 * Class:     nb_geometricBot_Model
 * Method:    substract
 * Signature: (Lnb/geometricBot/Model;Lnb/geometricBot/Model;)Lnb/geometricBot/Model;
 */
JNIEXPORT jobject JNICALL
Java_nb_geometricBot_Model_substract(JNIEnv *env, jclass class,
				     jobject jModelA, jobject jModelB);

/*
 * Class:     nb_geometricBot_Model
 * Method:    verify
 * Signature: ()Lnb/geometricBot/ModelStatus;
 */
JNIEXPORT jobject JNICALL
Java_nb_geometricBot_Model_verify(JNIEnv *env, jobject jModel);

#endif
