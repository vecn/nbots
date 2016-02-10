#ifndef __JNI_NB_GEOMETRIC_BOT_JMODEL_H__
#define __JNI_NB_GEOMETRIC_BOT_JMODEL_H__

#include <jni.h>

#include "nb/geometric_bot/model/model2D.h"

void get_model_from_java(JNIEnv *env, vcn_model_t *model, jobject jmodel);

jint jModel_getNVertices(JNIEnv *env, jobject jmodel);
jint jModel_getNEdges(JNIEnv *env, jobject jmodel);
jint jModel_getNHoles(JNIEnv *env, jobject jmodel);

#endif
