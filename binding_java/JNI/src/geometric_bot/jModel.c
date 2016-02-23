#include <stdlib.h>
#include <alloca.h>

#include <jni.h>

#include "nb/geometric_bot/model/model2D.h"

#include "geometric_bot/jModel.h"

#include "../src/nb/geometric_bot/model/model2D_struct.h"

static void create_vertices_from_model(JNIEnv *env, jclass class,
				       jobject jModel,
				       const vcn_model_t *model);
static void create_segments_from_model(JNIEnv *env, jclass class,
				       jobject jModel,
				       const vcn_model_t *model);
static void create_holes_from_model(JNIEnv *env, jclass class,
				       jobject jModel,
				       const vcn_model_t *model);
static void set_vertices_from_jModel(JNIEnv *env, float *vertices,
				     jint N, jobject jmodel);
static jfloatArray jModel_getVerticesRef(JNIEnv *env, jobject jmodel);
static void set_segments_from_jModel(JNIEnv *env, uint32_t *segments,
				     jint N, jobject jmodel);
static jintArray jModel_getEdgesRef(JNIEnv *env, jobject jmodel);
static void set_holes_from_jModel(JNIEnv *env, float *holes,
				  jint N, jobject jmodel);
static jfloatArray jModel_getHolesRef(JNIEnv *env, jobject jmodel);

jobject jModel_create(JNIEnv *env)
{
	const char *str_class = "nb/geometricBot/Model";
	jclass class = (*env)->FindClass(env, str_class);
	jmethodID method_id = (*env)->GetMethodID(env, class, "<init>", "()V");
	jobject instance = (*env)->NewObject(env, class, method_id);
	return instance;
}

void vcn_model_load_from_jModel(JNIEnv *env, vcn_model_t *model,
				const jobject jModel)
{
	jint N_vtx = jModel_getNVertices(env, jModel);
	jint N_sgm = jModel_getNEdges(env, jModel);
	jint N_holes = jModel_getNHoles(env, jModel);

	float *vertices = alloca(2 * N_vtx * sizeof(*vertices));
	set_vertices_from_jModel(env, vertices, N_vtx, jModel);

	uint32_t *segments = alloca(2 * N_sgm * sizeof(*segments));
	set_segments_from_jModel(env, segments, N_sgm, jModel);

	float *holes = alloca(2 * N_holes * sizeof(*holes));
	set_holes_from_jModel(env, holes, N_holes, jModel);

	vcn_model_load_from_farrays(model, vertices, N_vtx,
				    segments, N_sgm, holes, N_holes);
}

void jModel_load_from_model(JNIEnv *env, jobject jModel,
			      const vcn_model_t *model)
{
	jclass class =
		(*env)->GetObjectClass(env, jModel);
	create_vertices_from_model(env, class, jModel, model);
	create_segments_from_model(env, class, jModel, model);
	create_holes_from_model(env, class, jModel, model);
}

static void create_vertices_from_model(JNIEnv *env, jclass class,
				       jobject jModel,
				       const vcn_model_t *model)
{
	jfieldID field_id = 
		(*env)->GetFieldID(env, class, "vertices", "[F");
	jfloatArray jarray =
		(*env)->NewFloatArray(env, 2 * model->N);
	jfloat *array =
		(*env)->GetFloatArrayElements(env, jarray, NULL);

	for (uint32_t i = 0; i < 2 * model->N; i++)
		array[i] = model->vertex[i];
	(*env)->ReleaseFloatArrayElements(env, jarray, array, 0);
	(*env)->SetObjectField(env, jModel, field_id, jarray);
	
}

static void create_segments_from_model(JNIEnv *env, jclass class,
				       jobject jModel,
				       const vcn_model_t *model)
{
	if (0 < model->M) {
		jfieldID field_id = 
			(*env)->GetFieldID(env, class, "edges", "[I");
		jintArray jarray =
			(*env)->NewIntArray(env, 2 * model->M);
		jint *array =
			(*env)->GetIntArrayElements(env, jarray, NULL);

		for (uint32_t i = 0; i < 2 * model->M; i++)
			array[i] = model->edge[i];
		(*env)->ReleaseIntArrayElements(env, jarray, array, 0);
		(*env)->SetObjectField(env, jModel, field_id, jarray);
	}
}

static void create_holes_from_model(JNIEnv *env, jclass class,
				    jobject jModel,
				    const vcn_model_t *model)
{
	if (0 < model->H) {
		jfieldID field_id = 
			(*env)->GetFieldID(env, class, "holes", "[F");
		jfloatArray jarray =
			(*env)->NewFloatArray(env, 2 * model->H);
		jfloat *array =
			(*env)->GetFloatArrayElements(env, jarray, NULL);

		for (uint32_t i = 0; i < 2 * model->H; i++)
			array[i] = model->holes[i];
		(*env)->ReleaseFloatArrayElements(env, jarray, array, 0);
		(*env)->SetObjectField(env, jModel, field_id, jarray);
	}
}

jobject jModel_get_combination(JNIEnv *env, jclass class,
			       jobject jModelA, jobject jModelB,
			       void (*combine)(vcn_model_t *model,
					       const vcn_model_t *model1,
					       const vcn_model_t *model2,
					       double))
{
	vcn_model_t *model1 = alloca(vcn_model_get_memsize());
	vcn_model_init(model1);
	vcn_model_load_from_jModel(env, model1, jModelA);

	vcn_model_t *model2 = alloca(vcn_model_get_memsize());
	vcn_model_init(model2);
	vcn_model_load_from_jModel(env, model2, jModelB);

	vcn_model_t *model = alloca(vcn_model_get_memsize());
	vcn_model_init(model);
	combine(model, model1, model2, 0);

	jobject jModel = jModel_create(env);
	jModel_load_from_model(env, jModel, model);

	vcn_model_finish(model1);
	vcn_model_finish(model2);
	vcn_model_finish(model);

	return jModel;
}

jint jModel_getNVertices(JNIEnv *env, jobject jModel)
{
	jclass class = (*env)->GetObjectClass(env, jModel);
	jmethodID method_id =
		(*env)->GetMethodID(env, class, "getNVertices", "()I");
	jint N = (*env)->CallIntMethod(env, jModel, method_id);
	return N;
}

jint jModel_getNEdges(JNIEnv *env, jobject jModel)
{
	jclass class = (*env)->GetObjectClass(env, jModel);
	jmethodID method_id =
		(*env)->GetMethodID(env, class, "getNEdges", "()I");
	jint N = (*env)->CallIntMethod(env, jModel, method_id);
	return N;
}

jint jModel_getNHoles(JNIEnv *env, jobject jModel)
{
	jclass class = (*env)->GetObjectClass(env, jModel);
	jmethodID method_id = 
		(*env)->GetMethodID(env, class, "getNHoles", "()I");
	jint N = (*env)->CallIntMethod(env, jModel, method_id);
	return N;
}

static void set_vertices_from_jModel(JNIEnv *env, float *vertices, jint N,
				     jobject jModel)
{
	jfloatArray jvertices = jModel_getVerticesRef(env, jModel);
	jfloat *fvertices =
		(*env)->GetFloatArrayElements(env, jvertices, NULL);
	for (int i = 0; i < 2 * N; i++)
		vertices[i] = fvertices[i];
	(*env)->ReleaseFloatArrayElements(env, jvertices, fvertices, JNI_ABORT);
}

static jfloatArray jModel_getVerticesRef(JNIEnv *env, jobject jModel)
{
	jclass class = (*env)->GetObjectClass(env, jModel);
	jmethodID method_id =
		(*env)->GetMethodID(env, class, "getVerticesRef", "()[F");
	jfloatArray ref = (*env)->CallObjectMethod(env, jModel, method_id);
	return ref;	
}

static void set_segments_from_jModel(JNIEnv *env, uint32_t *segments, jint N,
				     jobject jModel)
{
	jintArray jsegments = jModel_getEdgesRef(env, jModel);
	jint *jedges =
		(*env)->GetIntArrayElements(env, jsegments, NULL);
	for (int i = 0; i < 2 * N; i++)
		segments[i] = jedges[i];
	(*env)->ReleaseIntArrayElements(env, jsegments, jedges, JNI_ABORT);
}

static jintArray jModel_getEdgesRef(JNIEnv *env, jobject jModel)
{
	jclass class = (*env)->GetObjectClass(env, jModel);
	jmethodID method_id =
		(*env)->GetMethodID(env, class, "getEdgesRef", "()[I");
	jintArray ref = (*env)->CallObjectMethod(env, jModel, method_id);
	return ref;	
}

static void set_holes_from_jModel(JNIEnv *env, float *holes, jint N,
				  jobject jModel)
{
	jfloatArray jholes = jModel_getHolesRef(env, jModel);
	jfloat *fholes =
		(*env)->GetFloatArrayElements(env, jholes, NULL);
	for (int i = 0; i < 2 * N; i++)
		holes[i] = fholes[i];
	(*env)->ReleaseFloatArrayElements(env, jholes, fholes, JNI_ABORT);
}

static jfloatArray jModel_getHolesRef(JNIEnv *env, jobject jModel)
{
	jclass class = (*env)->GetObjectClass(env, jModel);
	jmethodID method_id =
		(*env)->GetMethodID(env, class, "getHolesRef", "()[F");
	jfloatArray ref = (*env)->CallObjectMethod(env, jModel, method_id);
	return ref;	
}
