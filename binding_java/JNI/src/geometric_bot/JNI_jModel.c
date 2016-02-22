#include "geometric_bot/jModel.h"
#include "geometric_bot/JNI_jGeometricBot.h"

JNIEXPORT jobject JNICALL
Java_nb_geometricBot_Model_combine(JNIEnv *env, jclass class,
				   jobject jModelA, jobject jModelB)
{
	return NULL;
}

JNIEXPORT jobject JNICALL
Java_nb_geometricBot_Model_unify(JNIEnv *env, jclass class,
				 jobject jModelA, jobject jModelB)
{
	return NULL;
}

/*
 * Class:     nb_geometricBot_Model
 * Method:    intersect
 * Signature: (Lnb/geometricBot/Model;Lnb/geometricBot/Model;)Lnb/geometricBot/Model;
 */
JNIEXPORT jobject JNICALL
Java_nb_geometricBot_Model_intersect(JNIEnv *env, jclass class,
				     jobject jModelA, jobject jModelB)
{
	return NULL;
}

/*
 * Class:     nb_geometricBot_Model
 * Method:    substract
 * Signature: (Lnb/geometricBot/Model;Lnb/geometricBot/Model;)Lnb/geometricBot/Model;
 */
JNIEXPORT jobject JNICALL
Java_nb_geometricBot_Model_substract(JNIEnv *env, jclass class,
				     jobject jModelA, jobject jModelB)
{
	return NULL;
}

/*
 * Class:     nb_geometricBot_Model
 * Method:    verify
 * Signature: ()Lnb/geometricBot/ModelStatus;
 */
JNIEXPORT jobject JNICALL
Java_nb_geometricBot_Model_verify(JNIEnv *env, jobject jModel)
{
	return NULL;
}
