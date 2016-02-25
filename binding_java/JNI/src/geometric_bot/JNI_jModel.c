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

JNIEXPORT jobject JNICALL
Java_nb_geometricBot_Model_intersect(JNIEnv *env, jclass class,
				     jobject jModelA, jobject jModelB)
{
	return jModel_get_combination(env, class, jModelA, jModelB,
				      vcn_model_get_intersection);
}

JNIEXPORT jobject JNICALL
Java_nb_geometricBot_Model_substract(JNIEnv *env, jclass class,
				     jobject jModelA, jobject jModelB)
{
	return jModel_get_combination(env, class, jModelA, jModelB,
				      vcn_model_get_substraction);
}

JNIEXPORT jobject JNICALL
Java_nb_geometricBot_Model_difference(JNIEnv *env, jclass class,
				      jobject jModelA, jobject jModelB)
{
	return jModel_get_combination(env, class, jModelA, jModelB,
				      vcn_model_get_difference);
}

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

JNIEXPORT jboolean JNICALL
Java_nb_geometricBot_Model_isPointInside(JNIEnv *env, jobject jModel,
					 jdouble x, jdouble y)
{
	double vtx[2] = {x, y};

	vcn_model_t *model = alloca(vcn_model_get_memsize());
	vcn_model_init(model);
	vcn_model_load_from_jModel(env, model, jModel);

	jboolean is_inside = vcn_model_is_vtx_inside(model, vtx);

	vcn_model_finish(model);
	return is_inside;
}
