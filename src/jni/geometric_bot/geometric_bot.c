/******************************************************************************
 *   JNI Geometric Bot: Java Native Interface of the Geometric Bot            *
 *   2013-2015 Victor Eduardo Cardoso Nungaray                                *
 *   Twitter: @victore_cardoso                                                *
 *   email: victorc@cimat.mx                                                  *
 ******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "vcn/image_bot.h"
#include "vcn/geometric_bot.h"
#include "geometric_bot.h"

#include "../../vcn/geometric_bot/model/model2D_struct.h"

static vcn_model_t* create_model_from_JNI(JNIEnv *env,
					  jint N, jfloatArray jvertices,
					  jint M, jintArray jedges,
					  jint H, jfloatArray jholes);

static void model_load_from_jtypes(vcn_model_t *model,
				   jfloat vertices[], jint N_vertices,
				   jint segments[], jint N_segments,
				   jfloat holes[], jint N_holes);
static int load_mesh_into_JNI
               (JNIEnv *env, jclass class,
		const vcn_msh3trg_t *const msh3trg,
		jboolean include_edges, jboolean include_elems,
		jboolean include_elems_adj, 
		jboolean include_input_vtx, jboolean include_input_sgm,
		jobjectArray mesh_vtx, jobjectArray mesh_conn_edge,
		jobjectArray mesh_conn_mtx, jobjectArray mesh_conn_adj,
		jobjectArray mesh_input_vtx, jobjectArray mesh_input_sgm);

static int load_vertices_into_JNI(JNIEnv *env,
				  const vcn_msh3trg_t *const msh3trg,
				  jobjectArray jmesh_vtx);

static int load_edges_into_JNI(JNIEnv *env,
			       const vcn_msh3trg_t *const msh3trg,
			       jobjectArray jmesh_conn_edge);

static int load_connectivity_mtx_into_JNI(JNIEnv *env,
					  const vcn_msh3trg_t *const msh3trg,
					  jobjectArray jmesh_conn_mtx);

static int load_connectivity_adj_into_JNI(JNIEnv *env,
					  const vcn_msh3trg_t *const msh3trg,
					  jobjectArray jmesh_conn_adj);

static int load_input_vtx_into_JNI(JNIEnv *env,
				   const vcn_msh3trg_t *const msh3trg,
				   jobjectArray jmesh_input_vtx);

static int load_input_sgm_into_JNI(JNIEnv *env,
				   const vcn_msh3trg_t *const msh3trg,
				   jobjectArray jmesh_input_sgm);

JNIEXPORT jint JNICALL Java_mx_cimat_vcn_geometricBot_GeometricBot_JNICreateMesh__I_3FI_3II_3FIIFFFZZZZZ_3_3F_3_3I_3_3I_3_3I_3_3I_3_3_3I
                      (JNIEnv *env, jclass class,
		       jint N, jfloatArray vertices,
		       jint M, jintArray edges,
		       jint H, jfloatArray holes,
		       jint max_vtx, jint max_trg, jfloat min_angle,
		       jfloat max_edge_length, jfloat max_sgm_length,
		       jboolean include_edges, jboolean include_elems,
		       jboolean include_elems_adj, 
		       jboolean include_input_vtx, jboolean include_input_sgm,
		       jobjectArray mesh_vtx, jobjectArray mesh_conn_edge,
		       jobjectArray mesh_conn_mtx, jobjectArray mesh_conn_adj,
		       jobjectArray mesh_input_vtx, jobjectArray mesh_input_sgm)
{
	vcn_model_t *model =
		create_model_from_JNI(env, N, vertices, M, edges, H, holes);

	vcn_mesh_t *mesh = vcn_mesh_create();
	vcn_mesh_set_size_constraint(mesh, VCN_MESH_SIZE_CONSTRAINT_MAX_VTX,
				     max_vtx);
	vcn_mesh_set_size_constraint(mesh, VCN_MESH_SIZE_CONSTRAINT_MAX_TRG,
				     max_trg);
	vcn_mesh_set_geometric_constraint(mesh,
					  VCN_MESH_GEOM_CONSTRAINT_MIN_ANGLE,
					  min_angle);
	vcn_mesh_set_geometric_constraint(mesh,
					  VCN_MESH_GEOM_CONSTRAINT_MAX_EDGE_LENGTH,
					  max_edge_length);
	vcn_mesh_set_geometric_constraint(mesh,
					  VCN_MESH_GEOM_CONSTRAINT_MAX_SUBSGM_LENGTH,
					  max_sgm_length);
	vcn_mesh_generate_from_model(mesh, model);
	vcn_model_destroy(model);

	vcn_msh3trg_t *msh3trg =
		vcn_mesh_get_msh3trg(mesh,
				     include_edges,
				     include_elems,
				     include_elems_adj,
				     include_input_vtx,
				     include_input_sgm);
	vcn_mesh_destroy(mesh);

	if (NULL == msh3trg)
		return 2;

	int status = load_mesh_into_JNI(env, class, msh3trg, include_edges,
					include_elems, include_elems_adj,
					include_input_vtx, include_input_sgm,
					mesh_vtx, mesh_conn_edge, mesh_conn_mtx,
					mesh_conn_adj, mesh_input_vtx,
					mesh_input_sgm);
	if (0 != status)
		return 3;

	vcn_msh3trg_destroy(msh3trg);
	return 0;
}

JNIEXPORT jint JNICALL Java_mx_cimat_vcn_geometricBot_GeometricBot_JNICreateMesh__I_3FI_3II_3FLjava_lang_String_2IIFFFFZZZZZ_3_3F_3_3I_3_3I_3_3I_3_3I_3_3_3I
                      (JNIEnv *env, jclass class,
		       jint N, jfloatArray vertices,
		       jint M, jintArray edges,
		       jint H, jfloatArray holes,
		       jstring jimage_file, jint max_vtx, jint max_trg,
		       jfloat min_angle,
		       jfloat scale, jfloat x_disp, jfloat y_disp,
		       jboolean include_edges, jboolean include_elems,
		       jboolean include_elems_adj, 
		       jboolean include_input_vtx, jboolean include_input_sgm,
		       jobjectArray mesh_vtx, jobjectArray mesh_conn_edge, 
		       jobjectArray mesh_conn_mtx, jobjectArray mesh_conn_adj,
		       jobjectArray mesh_input_vtx, jobjectArray mesh_input_sgm)
/* (N = 0) indicates to that the model must be adjusted to the picture.
 */
{
	const char *image_file = (*env)->GetStringUTFChars(env, jimage_file, NULL);

	if(NULL == image_file)
		return 1;

	vcn_image_t *img = vcn_image_create();
	vcn_image_read(img, image_file);

	vcn_model_t *model;
	if (0 == N) {
		uint32_t w = vcn_image_get_width(img);
		uint32_t h = vcn_image_get_height(img);
		model = vcn_model_create_rectangle(0.0, 0.0, w, h);
	} else {
		model = create_model_from_JNI(env, N, vertices, M, edges, H, holes);
	}

	(*env)->ReleaseStringUTFChars(env, jimage_file, image_file);

	vcn_mesh_t *mesh = vcn_mesh_create();
	vcn_mesh_set_size_constraint(mesh, VCN_MESH_SIZE_CONSTRAINT_MAX_VTX,
				     max_vtx);
	vcn_mesh_set_size_constraint(mesh, VCN_MESH_SIZE_CONSTRAINT_MAX_TRG,
				     max_trg);
	vcn_mesh_set_geometric_constraint(mesh,
					  VCN_MESH_GEOM_CONSTRAINT_MIN_ANGLE,
					  min_angle);
	vcn_mesh_set_img_density(mesh, img, 0.0);
	vcn_mesh_generate_from_model(mesh, model);
	vcn_mesh_clear_img_density(mesh);
	vcn_image_destroy(img);
	vcn_model_destroy(model);

	vcn_msh3trg_t *msh3trg =
		vcn_mesh_get_msh3trg(mesh, include_edges, 
				     include_elems,
				     include_elems_adj, 
				     include_input_vtx,
				     include_input_sgm);

	vcn_mesh_destroy(mesh);
	if (NULL == msh3trg)
		return 3;

	int status = load_mesh_into_JNI(env, class, msh3trg, include_edges,
					include_elems, include_elems_adj,
					include_input_vtx, include_input_sgm,
					mesh_vtx, mesh_conn_edge, mesh_conn_mtx,
					mesh_conn_adj, mesh_input_vtx,
					mesh_input_sgm);
	if (0 != status)
		return 1;

	vcn_msh3trg_destroy(msh3trg);
	return 0;
}

static vcn_model_t* create_model_from_JNI(JNIEnv *env,
					  jint N, jfloatArray jvertices,
					  jint M, jintArray jedges,
					  jint H, jfloatArray jholes)
{
	jfloat *vertices = 
		(*env)->GetFloatArrayElements(env, jvertices, NULL);

	jint *edges = (*env)->GetIntArrayElements(env, jedges, NULL);

	jfloat *holes = NULL;
	if (0 < H)
		holes = (*env)->GetFloatArrayElements(env, jholes, NULL);
	vcn_model_t *model = vcn_model_create();
	model_load_from_jtypes(model, vertices, N, edges, M, holes, H);

	(*env)->ReleaseFloatArrayElements(env, jvertices, vertices, JNI_ABORT);
	(*env)->ReleaseIntArrayElements(env, jedges, edges, JNI_ABORT);
	if (0 < H)
		(*env)->ReleaseFloatArrayElements(env, jholes, holes, JNI_ABORT);

	return model;
}

static void model_load_from_jtypes(vcn_model_t *model,
				   jfloat vertices[], jint N_vertices,
				   jint segments[], jint N_segments,
				   jfloat holes[], jint N_holes)
{
	model->N = N_vertices;
	model->M = N_segments;
	model->H = N_holes;
	if (0 < N_vertices) {
		model_alloc_vertices(model);
		for (int i = 0; i < 2 * N_vertices; i++)
			model->vertex[i] = vertices[i];
	}
	if (0 < N_segments) {
		model_alloc_edges(model);
		for (int i = 0; i < 2 * N_segments; i++)
			model->edge[i] = segments[i];
	}
	if (0 < N_holes ) {
		model_alloc_holes(model);
		for (uint32_t i = 0; i < 2 * N_holes; i++)
			model->holes[i] = holes[i];
	}	
}

static int load_mesh_into_JNI
               (JNIEnv *env, jclass class,
		const vcn_msh3trg_t *const restrict msh3trg,
		jboolean include_edges, jboolean include_elems,
		jboolean include_elems_adj, 
		jboolean include_input_vtx, jboolean include_input_sgm,
		jobjectArray mesh_vtx, jobjectArray mesh_conn_edge,
		jobjectArray mesh_conn_mtx, jobjectArray mesh_conn_adj,
		jobjectArray mesh_input_vtx, jobjectArray mesh_input_sgm)
{
	int status = load_vertices_into_JNI(env, msh3trg, mesh_vtx);

	if (0 != status)
		return 1;

	if (include_edges) {
		status = load_edges_into_JNI(env, msh3trg, mesh_conn_edge);
		if (0 != status)
			return 2;
	}

	if (include_elems) {
		status = load_connectivity_mtx_into_JNI(env, msh3trg,
							mesh_conn_mtx);
		if (0 != status)
			return 3;
	}

	if (include_elems_adj) {
		status = load_connectivity_adj_into_JNI(env, msh3trg,
							mesh_conn_adj);
		if (0 != status)
			return 4;
	}

	if (include_input_vtx) {
		status = load_input_vtx_into_JNI(env, msh3trg, mesh_input_vtx);
		if (0 != status)
			return 5;
	}

	if (include_input_sgm) {
		status = load_input_sgm_into_JNI(env, msh3trg, mesh_input_sgm);
		if (0 != status)
			return 6;
	}
	return 0;
}

static int load_vertices_into_JNI(JNIEnv *env,
				  const vcn_msh3trg_t *const restrict msh3trg,
				  jobjectArray jmesh_vtx)
{
	jfloatArray jvertices =
		(*env)->NewFloatArray(env, 2 * msh3trg->N_vertices);

	if (NULL == jvertices) 
		return 1;

	jfloat *vertices =
		(*env)->GetFloatArrayElements(env, jvertices, NULL);

	for (uint32_t i = 0; i < 2 * msh3trg->N_vertices; i++)
		vertices[i] = (jfloat) msh3trg->vertices[i];

	(*env)->ReleaseFloatArrayElements(env, jvertices, vertices, 0);

	(*env)->SetObjectArrayElement(env, jmesh_vtx, 0, jvertices);
	return 0;
}

static int load_edges_into_JNI(JNIEnv *env,
			       const vcn_msh3trg_t *const restrict msh3trg,
			       jobjectArray jmesh_conn_edge)
{
	jintArray jconn_edge = (*env)->NewIntArray(env, msh3trg->N_edges * 2);

	if (NULL == jconn_edge) 
		return 1;

	jint *conn_edge = (*env)->GetIntArrayElements(env, jconn_edge, NULL);

	for (uint32_t i = 0; i < msh3trg->N_edges; i++) {
		conn_edge[i * 2] = (jint) msh3trg->edges[i * 2];
		conn_edge[i*2+1] = (jint) msh3trg->edges[i*2+1];
	}

	(*env)->ReleaseIntArrayElements(env, jconn_edge, conn_edge, 0);

	(*env)->SetObjectArrayElement(env, jmesh_conn_edge, 0, jconn_edge);  
	return 0;
}

static int load_connectivity_mtx_into_JNI
			(JNIEnv *env,
			 const vcn_msh3trg_t *const restrict msh3trg,
			 jobjectArray jmesh_conn_mtx)
{
	jintArray jconn_mtx = (*env)->NewIntArray(env, msh3trg->N_triangles * 3);

	if (NULL == jconn_mtx) 
		return 1;

	jint *conn_mtx = (*env)->GetIntArrayElements(env, jconn_mtx, NULL);

	for (uint32_t i = 0; i < msh3trg->N_triangles; i++) {
		for (uint32_t j = 0; j < 3; j++)
			conn_mtx[i*3+j] = (jint)
				msh3trg->vertices_forming_triangles[i*3+j];
	}
	(*env)->ReleaseIntArrayElements(env, jconn_mtx, conn_mtx, 0);
	(*env)->SetObjectArrayElement(env, jmesh_conn_mtx, 0, jconn_mtx);  
	return 0;
}


static int load_connectivity_adj_into_JNI(JNIEnv *env,
					  const vcn_msh3trg_t *const msh3trg,
					  jobjectArray jmesh_conn_adj)
{
	jintArray jconn_adj = (*env)->NewIntArray(env, msh3trg->N_triangles * 3);

	if (NULL == jconn_adj) 
		return 1;

	jint *conn_adj = (*env)->GetIntArrayElements(env, jconn_adj, NULL);

	for (uint32_t i = 0; i < msh3trg->N_triangles; i++) {
		for (uint32_t j = 0; j < 3; j++)
			conn_adj[i*3+j] = (jint)
				msh3trg->triangles_sharing_sides[i*3+j];
	}
	(*env)->ReleaseIntArrayElements(env, jconn_adj, conn_adj, 0);
	(*env)->SetObjectArrayElement(env, jmesh_conn_adj, 0, jconn_adj);
	return 0;
}

static int load_input_vtx_into_JNI(JNIEnv *env,
				   const vcn_msh3trg_t *const msh3trg,
				   jobjectArray jmesh_input_vtx)
{
	jintArray jinput_vtx = (*env)->NewIntArray(env, msh3trg->N_input_vertices);
	if (NULL == jinput_vtx) 
		return 1;
	
	jint *input_vtx = (*env)->GetIntArrayElements(env, jinput_vtx, NULL);

	for (uint32_t i = 0; i < msh3trg->N_input_vertices; i++)
		input_vtx[i] = (jint) msh3trg->input_vertices[i];

	(*env)->ReleaseIntArrayElements(env, jinput_vtx, input_vtx, 0);
	(*env)->SetObjectArrayElement(env, jmesh_input_vtx, 0, jinput_vtx);  
	return 0;
}

static int load_input_sgm_into_JNI(JNIEnv *env,
				   const vcn_msh3trg_t *const msh3trg,
				   jobjectArray jmesh_input_sgm)
{
	jclass class = (*env)->FindClass(env, "[I");
	jobjectArray jarray_input_sgm = 
		(*env)->NewObjectArray(env, msh3trg->N_input_segments, class, 0);

	if (NULL == jarray_input_sgm) 
		return 1;

	for (uint32_t i = 0; i < msh3trg->N_input_segments; i++) {
		jintArray jinput_sgm =
			(*env)->NewIntArray(env,
					    msh3trg->N_subsgm_x_inputsgm[i] + 1);
		if (NULL == jinput_sgm) 
			return 2;

		jint *input_sgm = (*env)->GetIntArrayElements(env, 
							      jinput_sgm, NULL);

		for (uint32_t j = 0; j < msh3trg->N_subsgm_x_inputsgm[i] + 1; j++)
			input_sgm[j] = (jint) msh3trg->meshvtx_x_inputsgm[i][j];

		(*env)->ReleaseIntArrayElements(env, jinput_sgm, input_sgm, 0);
		(*env)->SetObjectArrayElement(env, jarray_input_sgm, 0, jinput_sgm);
	}
	(*env)->SetObjectArrayElement(env, jmesh_input_sgm, 0, jarray_input_sgm);
	return 0;
}
