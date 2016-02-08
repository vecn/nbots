/**
 * This file is directly associated with the class
 *   vcn.pdeBot.finiteElement.solidMechanics.StaticElasticity2D
 *
 * This file was generated using javah:
 *   javah vcn.pdeBot.finiteElement.solidMechanics.StaticElasticity2D
 */

#ifndef __JNI_VCN_PDE_BOT_H__
#define __JNI_VCN_PDE_BOT_H__

#include <jni.h>

/*
 * Class:     vcn_pdeBot_finiteElement_solidMechanics_StaticElasticity2D
 * Method:    solve
 * Signature: (Lvcn/geometricBot/Model;Lvcn/pdeBot/Material;Lvcn/pdeBot/BoundaryConditions;Lvcn/pdeBot/BoundaryConditions;D)Lvcn/pdeBot/finiteElement/MeshResults;
 */
JNIEXPORT jobject JNICALL Java_vcn_pdeBot_finiteElement_solidMechanics_StaticElasticity2D_solve
		(JNIEnv *env, jclass class, jobject model, jobject material,
		 jobject bCondVtx, jobject bCondSgm, jdouble thickness);

#endif
