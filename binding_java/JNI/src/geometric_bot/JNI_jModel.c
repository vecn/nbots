#include <alloca.h>

#include <jni.h>

#include "nb/geometric_bot/model/model2D.h"
#include "nb/geometric_bot/model/modules2D/verifier.h"
#include "nb/geometric_bot/model/modules2D/clipper.h"

#include "geometric_bot/jModel.h"
#include "geometric_bot/jModelStatus.h"
#include "geometric_bot/JNI_jGeometricBot.h"

JNIEXPORT jobject JNICALL
Java_nb_geometricBot_Model_combine(JNIEnv *env, jclass class,
				   jobject jModelA, jobject jModelB)
{
	return jModel_get_combination(env, class, jModelA, jModelB,
				      vcn_model_get_combination);
}

JNIEXPORT jobject JNICALL
Java_nb_geometricBot_Model_unify(JNIEnv *env, jclass class,
				 jobject jModelA, jobject jModelB)
{
	return jModel_get_combination(env, class, jModelA, jModelB,
				      vcn_model_get_union);
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
	return jModel_get_combination(env, class, jModelA, jModelB,
				      vcn_model_get_intersection);
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
	return jModel_get_combination(env, class, jModelA, jModelB,
				      vcn_model_get_substraction);
}

/*
 * Class:     nb_geometricBot_Model
 * Method:    difference
 * Signature: (Lnb/geometricBot/Model;Lnb/geometricBot/Model;)Lnb/geometricBot/Model;
 */
JNIEXPORT jobject JNICALL
Java_nb_geometricBot_Model_difference(JNIEnv *env, jclass class,
				      jobject jModelA, jobject jModelB)
{
	return jModel_get_combination(env, class, jModelA, jModelB,
				      vcn_model_get_difference);
}

/*
 * Class:     nb_geometricBot_Model
 * Method:    verify
 * Signature: ()Lnb/geometricBot/ModelStatus;
 */
JNIEXPORT jobject JNICALL
Java_nb_geometricBot_Model_verify(JNIEnv *env, jobject jModel)
{	
	vcn_model_t *model = alloca(vcn_model_get_memsize());
	vcn_model_init(model);
	vcn_model_load_from_jModel(env, model, jModel);

	int status = vcn_model_verify_consistence(model, NULL);

	vcn_model_finish(model);

	jobject jModelStatus = jModelStatus_create(env, status);
	return jModelStatus;
}
