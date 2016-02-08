/**
 * This file is directly associated with the class
 *   vcn.geometricBot.GeometricBot
 *
 * This file was generated using javah:
 *   javah vcn.geometricBot.GeometricBot
 */


#ifndef __JNI_VCN_GEOMETRIC_BOT_H__
#define __JNI_VCN_GEOMETRIC_BOT_H__

#include <jni.h>

/*
 * Class:     vcn_geometricBot_GeometricBot
 * Method:    JNICreateMesh
 * Signature: (I[FI[II[FIIFFFZZZZZ[[F[[I[[I[[I[[I[[[I)I
 */

JNIEXPORT jint JNICALL Java_vcn_geometricBot_GeometricBot_JNICreateMesh__I_3FI_3II_3FIIFFFZZZZZ_3_3F_3_3I_3_3I_3_3I_3_3I_3_3_3I
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
		 jobjectArray mesh_input_vtx, jobjectArray mesh_input_sgm);

/*
 * Class:     vcn_geometricBot_GeometricBot
 * Method:    JNICreateMesh
 * Signature: (I[FI[II[FLjava/lang/String;IIFFFFZZZZZ[[F[[I[[I[[I[[I[[[I)I
 */

JNIEXPORT jint JNICALL Java_vcn_geometricBot_GeometricBot_JNICreateMesh__I_3FI_3II_3FLjava_lang_String_2IIFFFFZZZZZ_3_3F_3_3I_3_3I_3_3I_3_3I_3_3_3I
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
		 jobjectArray mesh_input_vtx, jobjectArray mesh_input_sgm);

#endif
