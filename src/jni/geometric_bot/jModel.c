#include <stdlib.h>
#include <jni.h>

#include "vcn/geometric_bot/model/model2D.h"

#include "jModel.h"

static void set_vertices_from_jModel(JNIEnv *env, float *vertices,
				     jint N, jobject jmodel);
static jfloatArray jModel_getVerticesRef(JNIEnv *env, jobject jmodel);
static void set_segments_from_jModel(JNIEnv *env, uint32_t *segments,
				     jint N, jobject jmodel);
static jintArray jModel_getEdgesRef(JNIEnv *env, jobject jmodel);
static void set_holes_from_jModel(JNIEnv *env, float *holes,
				  jint N, jobject jmodel);
static jfloatArray jModel_getHolesRef(JNIEnv *env, jobject jmodel);

void get_model_from_java(JNIEnv *env, vcn_model_t *model, jobject jmodel)
{
	jint N_vtx = jModel_getNVertices(env, jmodel);
	jint N_sgm = jModel_getNEdges(env, jmodel);
	jint N_holes = jModel_getNHoles(env, jmodel);
	float *vertices = malloc(2 * N_vtx * sizeof(*vertices));
	set_vertices_from_jModel(env, vertices, N_vtx, jmodel);
	uint32_t *segments = malloc(2 * N_sgm * sizeof(*segments));
	set_segments_from_jModel(env, segments, N_sgm, jmodel);
	float *holes = malloc(2 * N_holes * sizeof(*holes));
	set_holes_from_jModel(env, holes, N_holes, jmodel);

	vcn_model_load_from_farrays(model, vertices, N_vtx,
				    segments, N_sgm, holes, N_holes);
	free(vertices);
	free(segments);
	free(holes);
}

jint jModel_getNVertices(JNIEnv *env, jobject jmodel)
{
	jclass class = (*env)->GetObjectClass(env, jmodel);
	jmethodID method_id =
		(*env)->GetMethodID(env, class, "getNVertices", "()I");
	jint N = (*env)->CallIntMethod(env, jmodel, method_id);
	return N;
}

jint jModel_getNEdges(JNIEnv *env, jobject jmodel)
{
	jclass class = (*env)->GetObjectClass(env, jmodel);
	jmethodID method_id =
		(*env)->GetMethodID(env, class, "getNEdges", "()I");
	jint N = (*env)->CallIntMethod(env, jmodel, method_id);
	return N;
}

jint jModel_getNHoles(JNIEnv *env, jobject jmodel)
{
	jclass class = (*env)->GetObjectClass(env, jmodel);
	jmethodID method_id = 
		(*env)->GetMethodID(env, class, "getNHoles", "()I");
	jint N = (*env)->CallIntMethod(env, jmodel, method_id);
	return N;
}

static void set_vertices_from_jModel(JNIEnv *env, float *vertices, jint N,
				     jobject jmodel)
{
	jfloatArray jvertices = jModel_getVerticesRef(env, jmodel);
	jfloat *fvertices =
		(*env)->GetFloatArrayElements(env, jvertices, NULL);
	for (int i = 0; i < 2 * N; i++)
		vertices[i] = fvertices[i];
	(*env)->ReleaseFloatArrayElements(env, jvertices, fvertices, JNI_ABORT);
}

static jfloatArray jModel_getVerticesRef(JNIEnv *env, jobject jmodel)
{
	jclass class = (*env)->GetObjectClass(env, jmodel);
	jmethodID method_id =
		(*env)->GetMethodID(env, class, "getVerticesRef", "()[F");
	jfloatArray ref = (*env)->CallObjectMethod(env, jmodel, method_id);
	return ref;	
}

static void set_segments_from_jModel(JNIEnv *env, uint32_t *segments, jint N,
				     jobject jmodel)
{
	jintArray jsegments = jModel_getEdgesRef(env, jmodel);
	jint *jedges =
		(*env)->GetIntArrayElements(env, jsegments, NULL);
	for (int i = 0; i < 2 * N; i++)
		segments[i] = jedges[i];
	(*env)->ReleaseIntArrayElements(env, jsegments, jedges, JNI_ABORT);
}

static jintArray jModel_getEdgesRef(JNIEnv *env, jobject jmodel)
{
	jclass class = (*env)->GetObjectClass(env, jmodel);
	jmethodID method_id =
		(*env)->GetMethodID(env, class, "getEdgesRef", "()[I");
	jintArray ref = (*env)->CallObjectMethod(env, jmodel, method_id);
	return ref;	
}

static void set_holes_from_jModel(JNIEnv *env, float *holes, jint N,
				  jobject jmodel)
{
	jfloatArray jholes = jModel_getHolesRef(env, jmodel);
	jfloat *fholes =
		(*env)->GetFloatArrayElements(env, jholes, NULL);
	for (int i = 0; i < 2 * N; i++)
		holes[i] = fholes[i];
	(*env)->ReleaseFloatArrayElements(env, jholes, fholes, JNI_ABORT);
}

static jfloatArray jModel_getHolesRef(JNIEnv *env, jobject jmodel)
{
	jclass class = (*env)->GetObjectClass(env, jmodel);
	jmethodID method_id =
		(*env)->GetMethodID(env, class, "getHolesRef", "()[F");
	jfloatArray ref = (*env)->CallObjectMethod(env, jmodel, method_id);
	return ref;	
}
