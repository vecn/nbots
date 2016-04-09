#ifndef __JNI_NB_GEOMETRIC_BOT_LOAD_JMESH_H__
#define __JNI_NB_GEOMETRIC_BOT_LOAD_JMESH_H__

#include <jni.h>

#include "nb/geometric_bot/mesh/elements2D/triangles.h"

jobject jMesh_new(JNIEnv *env);
void load_jMesh_from_msh3trg(JNIEnv *env, const vcn_msh3trg_t *const msh3trg,
			     jobject jmesh);

#endif
