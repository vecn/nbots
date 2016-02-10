/**
 * This file is directly associated with the class
 *   vcn.pdeBot.finiteElement.solidMechanics.StaticElasticity2D
 *
 * This file was generated using javah:
 *   javah vcn.pdeBot.finiteElement.solidMechanics.StaticElasticity2D
 */

#ifndef __JNI_NB_PDE_BOT_H__
#define __JNI_NB_PDE_BOT_H__

#include <jni.h>

/*
 * Class:     nb_pdeBot_finiteElement_solidMechanics_StaticElasticity2D
 * Method:    solve
 * Signature: (Lnb/geometricBot/Model;Lnb/pdeBot/Material;Lnb/pdeBot/BoundaryConditions;Lnb/pdeBot/BoundaryConditions;Lnb/pdeBot/BoundaryConditions;Lnb/pdeBot/BoundaryConditions;D)Lnb/pdeBot/finiteElement/MeshResults;
 */
JNIEXPORT jobject JNICALL
Java_nb_pdeBot_finiteElement_solidMechanics_StaticElasticity2D_solve
		(JNIEnv *env, jclass class, jobject jmodel, jobject jmaterial,
		 jobject jBCDirichletVtx, jobject jBCNeumannVtx,
		 jobject jBCDirichletSgm, jobject jBCNeumannSgm,
		 jdouble jthickness);

#endif
