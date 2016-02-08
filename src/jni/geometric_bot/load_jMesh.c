#include <jni.h>

#include "vcn/geometric_bot/mesh/elements2D/triangles.h"

static void load_vertices_into_JNI(JNIEnv *env,
				   const vcn_msh3trg_t *const restrict msh3trg,
				   jobjectArray jmesh_vtx);
static void load_edges_into_JNI(JNIEnv *env,
				const vcn_msh3trg_t *const restrict msh3trg,
				jobjectArray jmesh_conn_edge);
static void load_connectivity_mtx_into_JNI
			(JNIEnv *env,
			 const vcn_msh3trg_t *const restrict msh3trg,
			 jobjectArray jmesh_conn_mtx);
static void load_connectivity_adj_into_JNI(JNIEnv *env,
					   const vcn_msh3trg_t *const msh3trg,
					   jobjectArray jmesh_conn_adj);
static void load_input_vtx_into_JNI(JNIEnv *env,
				    const vcn_msh3trg_t *const msh3trg,
				    jobjectArray jmesh_input_vtx);
static void load_input_sgm_into_JNI(JNIEnv *env,
				    const vcn_msh3trg_t *const msh3trg,
				    jobjectArray jmesh_input_sgm);

void load_jMesh_from_msh3trg(JNIEnv *env, const vcn_msh3trg_t *const msh3trg,
			     jobject jmesh)
{
	jfloatArray jmesh_vtx;
	load_vertices_into_JNI(env, msh3trg, mesh_vtx);
	/* AQUI VOY */

	load_edges_into_JNI(env, msh3trg, mesh_conn_edge);
	load_connectivity_mtx_into_JNI(env, msh3trg, mesh_conn_mtx);
	load_connectivity_adj_into_JNI(env, msh3trg, mesh_conn_adj);
	load_input_vtx_into_JNI(env, msh3trg, mesh_input_vtx);
	load_input_sgm_into_JNI(env, msh3trg, mesh_input_sgm);	
}

static void load_vertices_into_JNI(JNIEnv *env,
				   const vcn_msh3trg_t *const restrict msh3trg,
				   jobjectArray jmesh_vtx)
{
	jfloatArray jvertices =
		(*env)->NewFloatArray(env, 2 * msh3trg->N_vertices);

	jfloat *vertices =
		(*env)->GetFloatArrayElements(env, jvertices, NULL);

	for (uint32_t i = 0; i < 2 * msh3trg->N_vertices; i++)
		vertices[i] = (jfloat) msh3trg->vertices[i];

	(*env)->ReleaseFloatArrayElements(env, jvertices, vertices, 0);

	(*env)->SetObjectArrayElement(env, jmesh_vtx, 0, jvertices);
}

static void load_edges_into_JNI(JNIEnv *env,
				const vcn_msh3trg_t *const restrict msh3trg,
				jobjectArray jmesh_conn_edge)
{
	jintArray jconn_edge = (*env)->NewIntArray(env, msh3trg->N_edges * 2);

	jint *conn_edge = (*env)->GetIntArrayElements(env, jconn_edge, NULL);

	for (uint32_t i = 0; i < msh3trg->N_edges; i++) {
		conn_edge[i * 2] = (jint) msh3trg->edges[i * 2];
		conn_edge[i*2+1] = (jint) msh3trg->edges[i*2+1];
	}

	(*env)->ReleaseIntArrayElements(env, jconn_edge, conn_edge, 0);

	(*env)->SetObjectArrayElement(env, jmesh_conn_edge, 0, jconn_edge);  
}

static void load_connectivity_mtx_into_JNI
			(JNIEnv *env,
			 const vcn_msh3trg_t *const restrict msh3trg,
			 jobjectArray jmesh_conn_mtx)
{
	jintArray jconn_mtx = (*env)->NewIntArray(env, msh3trg->N_triangles * 3);

	jint *conn_mtx = (*env)->GetIntArrayElements(env, jconn_mtx, NULL);

	for (uint32_t i = 0; i < msh3trg->N_triangles; i++) {
		for (uint32_t j = 0; j < 3; j++)
			conn_mtx[i*3+j] = (jint)
				msh3trg->vertices_forming_triangles[i*3+j];
	}
	(*env)->ReleaseIntArrayElements(env, jconn_mtx, conn_mtx, 0);
	(*env)->SetObjectArrayElement(env, jmesh_conn_mtx, 0, jconn_mtx);
}


static void load_connectivity_adj_into_JNI(JNIEnv *env,
					   const vcn_msh3trg_t *const msh3trg,
					   jobjectArray jmesh_conn_adj)
{
	jintArray jconn_adj = (*env)->NewIntArray(env, msh3trg->N_triangles * 3);

	jint *conn_adj = (*env)->GetIntArrayElements(env, jconn_adj, NULL);

	for (uint32_t i = 0; i < msh3trg->N_triangles; i++) {
		for (uint32_t j = 0; j < 3; j++)
			conn_adj[i*3+j] = (jint)
				msh3trg->triangles_sharing_sides[i*3+j];
	}
	(*env)->ReleaseIntArrayElements(env, jconn_adj, conn_adj, 0);
	(*env)->SetObjectArrayElement(env, jmesh_conn_adj, 0, jconn_adj);
}

static void load_input_vtx_into_JNI(JNIEnv *env,
				    const vcn_msh3trg_t *const msh3trg,
				    jobjectArray jmesh_input_vtx)
{
	jintArray jinput_vtx = (*env)->NewIntArray(env, msh3trg->N_input_vertices);
	
	jint *input_vtx = (*env)->GetIntArrayElements(env, jinput_vtx, NULL);

	for (uint32_t i = 0; i < msh3trg->N_input_vertices; i++)
		input_vtx[i] = (jint) msh3trg->input_vertices[i];

	(*env)->ReleaseIntArrayElements(env, jinput_vtx, input_vtx, 0);
	(*env)->SetObjectArrayElement(env, jmesh_input_vtx, 0, jinput_vtx);
}

static void load_input_sgm_into_JNI(JNIEnv *env,
				    const vcn_msh3trg_t *const msh3trg,
				    jobjectArray jmesh_input_sgm)
{
	jclass class = (*env)->FindClass(env, "[I");
	jobjectArray jarray_input_sgm = 
		(*env)->NewObjectArray(env, msh3trg->N_input_segments, class, 0);

	for (uint32_t i = 0; i < msh3trg->N_input_segments; i++) {
		jintArray jinput_sgm =
			(*env)->NewIntArray(env,
					    msh3trg->N_subsgm_x_inputsgm[i] + 1);

		jint *input_sgm = (*env)->GetIntArrayElements(env, 
							      jinput_sgm, NULL);

		for (uint32_t j = 0; j < msh3trg->N_subsgm_x_inputsgm[i] + 1; j++)
			input_sgm[j] = (jint) msh3trg->meshvtx_x_inputsgm[i][j];

		(*env)->ReleaseIntArrayElements(env, jinput_sgm, input_sgm, 0);
		(*env)->SetObjectArrayElement(env, jarray_input_sgm, 0, jinput_sgm);
	}
	(*env)->SetObjectArrayElement(env, jmesh_input_sgm, 0, jarray_input_sgm);
}
