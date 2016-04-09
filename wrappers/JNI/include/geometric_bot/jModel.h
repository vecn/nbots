#ifndef __JNI_NB_GEOMETRIC_BOT_JMODEL_H__
#define __JNI_NB_GEOMETRIC_BOT_JMODEL_H__

#include <jni.h>

#include "nb/geometric_bot/model/model2D.h"

void vcn_model_load_from_jModel(JNIEnv *env, vcn_model_t *model,
				const jobject jModel);

jobject jModel_create(JNIEnv *env);
void jModel_load_from_model(JNIEnv *env, jobject jModel,
			      const vcn_model_t *model);
jobject jModel_get_combination(JNIEnv *env, jclass class,
			       jobject jModelA, jobject jModelB,
			       void (*combine)(vcn_model_t *model,
					       const vcn_model_t *model1,
					       const vcn_model_t *model2,
					       double));

jint jModel_getNVertices(JNIEnv *env, jobject jModel);
jint jModel_getNEdges(JNIEnv *env, jobject jModel);
jint jModel_getNHoles(JNIEnv *env, jobject jModel);

#endif
