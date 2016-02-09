#include <jni.h>

#include "vcn/geometric_bot/mesh/elements2D/triangles.h"

static void load_vertices_into_JNI(JNIEnv *env,
				   const vcn_msh3trg_t *const restrict msh3trg,
				   jobject jmesh);
static void load_edges_into_JNI(JNIEnv *env,
				const vcn_msh3trg_t *const restrict msh3trg,
				jobject jmesh);
static void load_connectivity_mtx_into_JNI
			(JNIEnv *env,
			 const vcn_msh3trg_t *const restrict msh3trg,
			 jobject jmesh);
static void load_connectivity_adj_into_JNI(JNIEnv *env,
					   const vcn_msh3trg_t *const msh3trg,
					   jobject jmesh);
static void load_input_vtx_into_JNI(JNIEnv *env,
				    const vcn_msh3trg_t *const msh3trg,
				    jobject jmesh);
static void load_input_sgm_into_JNI(JNIEnv *env,
				    const vcn_msh3trg_t *const msh3trg,
				    jobject jmesh);

jobject jMesh_new(JNIEnv *env)
{
	jclass class = (*env)->FindClass(env, "Lvcn/geometricBot/Mesh;");
	jmethodID method_id = (*env)->GetMethodID(env, class, "<init>", "()V");
	jobject instance = (*env)->NewObject(env, class, method_id);
	return instance;

}

void load_jMesh_from_msh3trg(JNIEnv *env, const vcn_msh3trg_t *const msh3trg,
			     jobject jmesh)
{
	load_vertices_into_JNI(env, msh3trg, jmesh);
	load_edges_into_JNI(env, msh3trg, jmesh);
	load_connectivity_mtx_into_JNI(env, msh3trg, jmesh);
	load_connectivity_adj_into_JNI(env, msh3trg, jmesh);
	load_input_vtx_into_JNI(env, msh3trg, jmesh);
	load_input_sgm_into_JNI(env, msh3trg, jmesh);	
}

static void load_vertices_into_JNI(JNIEnv *env,
				   const vcn_msh3trg_t *const restrict msh3trg,
				   jobject jmesh)
{
	jclass class =
		(*env)->GetObjectClass(env, jmesh);
	jfieldID field_id = 
		(*env)->GetFieldID(env, class, "vertices", "[F");
	jfloatArray jvertices =
		(*env)->NewFloatArray(env, 2 * msh3trg->N_vertices);
	jfloat *vertices =
		(*env)->GetFloatArrayElements(env, jvertices, NULL);

	for (uint32_t i = 0; i < 2 * msh3trg->N_vertices; i++)
		vertices[i] = msh3trg->vertices[i];

	(*env)->ReleaseFloatArrayElements(env, jvertices, vertices, 0);

	(*env)->SetObjectField(env, jmesh, field_id, jvertices);
}

static void load_edges_into_JNI(JNIEnv *env,
				const vcn_msh3trg_t *const restrict msh3trg,
				jobject jmesh)
{
	jclass class =
		(*env)->GetObjectClass(env, jmesh);
	jfieldID field_id = 
		(*env)->GetFieldID(env, class, "connEdges", "[I");
	jintArray jedges =
		(*env)->NewIntArray(env, 2 * msh3trg->N_edges);
	jint *edges =
		(*env)->GetIntArrayElements(env, jedges, NULL);

	for (uint32_t i = 0; i < msh3trg->N_edges; i++) {
		edges[i * 2] = msh3trg->edges[i * 2];
		edges[i*2+1] = msh3trg->edges[i*2+1];
	}

	(*env)->ReleaseIntArrayElements(env, jedges, edges, 0);
	(*env)->SetObjectField(env, jmesh, field_id, jedges);
}

static void load_connectivity_mtx_into_JNI
			(JNIEnv *env,
			 const vcn_msh3trg_t *const restrict msh3trg,
			 jobject jmesh)
{
	jclass class =
		(*env)->GetObjectClass(env, jmesh);
	jfieldID field_id = 
		(*env)->GetFieldID(env, class, "connMtx", "[I");
	jintArray jconnMtx =
		(*env)->NewIntArray(env, msh3trg->N_triangles * 3);
	jint *connMtx =
		(*env)->GetIntArrayElements(env, jconnMtx, NULL);

	for (uint32_t i = 0; i < msh3trg->N_triangles; i++) {
		for (uint32_t j = 0; j < 3; j++)
			connMtx[i*3+j] =
				msh3trg->vertices_forming_triangles[i*3+j];
	}
	(*env)->ReleaseIntArrayElements(env, jconnMtx, connMtx, 0);
	(*env)->SetObjectField(env, jmesh, field_id, jconnMtx);
}


static void load_connectivity_adj_into_JNI(JNIEnv *env,
					   const vcn_msh3trg_t *const msh3trg,
					   jobject jmesh)
{
	jclass class =
		(*env)->GetObjectClass(env, jmesh);
	jfieldID field_id = 
		(*env)->GetFieldID(env, class, "connAdj", "[I");
	jintArray jconnAdj =
		(*env)->NewIntArray(env, msh3trg->N_triangles * 3);
	jint *connAdj = (*env)->GetIntArrayElements(env, jconnAdj, NULL);

	for (uint32_t i = 0; i < msh3trg->N_triangles; i++) {
		for (uint32_t j = 0; j < 3; j++)
			connAdj[i*3+j] =
				msh3trg->triangles_sharing_sides[i*3+j];
	}
	(*env)->ReleaseIntArrayElements(env, jconnAdj, connAdj, 0);
	(*env)->SetObjectField(env, jmesh, field_id, jconnAdj);
}

static void load_input_vtx_into_JNI(JNIEnv *env,
				    const vcn_msh3trg_t *const msh3trg,
				    jobject jmesh)
{
	jclass class =
		(*env)->GetObjectClass(env, jmesh);
	jfieldID field_id = 
		(*env)->GetFieldID(env, class, "modelVtx", "[I");
	jintArray jmodelVtx =
		(*env)->NewIntArray(env, msh3trg->N_input_vertices);
	jint *modelVtx =
		(*env)->GetIntArrayElements(env, jmodelVtx, NULL);

	for (uint32_t i = 0; i < msh3trg->N_input_vertices; i++)
		modelVtx[i] = msh3trg->input_vertices[i];

	(*env)->ReleaseIntArrayElements(env, jmodelVtx, modelVtx, 0);
	(*env)->SetObjectField(env, jmesh, field_id, jmodelVtx);
}

static void load_input_sgm_into_JNI(JNIEnv *env,
				    const vcn_msh3trg_t *const msh3trg,
				    jobject jmesh)
{
	jclass jmesh_class =
		(*env)->GetObjectClass(env, jmesh);
	jfieldID field_id = 
		(*env)->GetFieldID(env, jmesh_class, "modelSgm", "[[I");

	jclass class = (*env)->FindClass(env, "[I");
	jobjectArray jmodelSgm = 
		(*env)->NewObjectArray(env, msh3trg->N_input_segments, class, 0);

	for (uint32_t i = 0; i < msh3trg->N_input_segments; i++) {
		jint N_vtx = msh3trg->N_subsgm_x_inputsgm[i] + 1;
		jintArray jinput_sgm = (*env)->NewIntArray(env, N_vtx);
		jint *input_sgm = 
			(*env)->GetIntArrayElements(env, jinput_sgm, NULL);
		for (uint32_t j = 0; j < N_vtx; j++)
			input_sgm[j] = msh3trg->meshvtx_x_inputsgm[i][j];

		(*env)->ReleaseIntArrayElements(env, jinput_sgm, input_sgm, 0);
		(*env)->SetObjectArrayElement(env, jmodelSgm, i, jinput_sgm);
	}
	(*env)->SetObjectField(env, jmesh, field_id, jmodelSgm);
}
