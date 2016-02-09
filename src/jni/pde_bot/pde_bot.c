#include <stdlib.h>

#include "vcn/eigen_bot.h"
#include "vcn/geometric_bot.h"
#include "vcn/pde_bot.h"

#include "pde_bot.h"
#include "../geometric_bot/load_jMesh.h"
#include "../geometric_bot/jModel.h"

static void get_bconditions_from_java(JNIEnv *env, vcn_bcond_t *bcond,
				      jobject jbCondVtx, jobject jbCondSgm);

static jint jBoundaryConditions_getN(JNIEnv *env, jobject jBoundaryConditions);
static jintArray jBoundaryConditions_getIdsRef(JNIEnv *env,
					       jobject jBoundaryConditions);
static jcharArray jBoundaryConditions_getDofRef
			(JNIEnv *env, jobject jBoundaryConditions);
static jdoubleArray jBoundaryConditions_getValuesRef
			(JNIEnv *env, jobject jBoundaryConditions);
static void set_bc_dof(bool enabled[2], jchar dof);
static void get_material_from_java(JNIEnv *env, vcn_fem_material_t *material,
				   jobject jmaterial);
static jdouble jMaterial_getPoissonModulus(JNIEnv *env, jobject jmaterial);
static jdouble jMaterial_getYoungModulus(JNIEnv *env, jobject jmaterial);
static jobject jMeshResults_new(JNIEnv *env);
static void set_results_into_jmesh(JNIEnv *env, jobject jmesh,
				   const double *displacement,
				   const double *strain);

/*
 * Class:     vcn_pdeBot_finiteElement_solidMechanics_StaticElasticity2D
 * Method:    solve
 * Signature: (Lvcn/geometricBot/Model;Lvcn/pdeBot/Material;Lvcn/pdeBot/BoundaryConditions;Lvcn/pdeBot/BoundaryConditions;D)Lvcn/pdeBot/finiteElement/MeshResults;
 */
JNIEXPORT jobject JNICALL
Java_vcn_pdeBot_finiteElement_solidMechanics_StaticElasticity2D_solve
		(JNIEnv *env, jclass class, jobject jmodel, jobject jmaterial,
		 jobject jbCondVtx, jobject jbCondSgm, jdouble jthickness)
{
	/* Read optimization parameters */
	vcn_bcond_t* bconditions = vcn_fem_bcond_create();
	get_bconditions_from_java(env, bconditions, jbCondVtx, jbCondSgm);
	
	vcn_fem_material_t* material = vcn_fem_material_create();
	get_material_from_java(env, material, jmaterial);
	
	vcn_model_t* model = vcn_model_create();
	get_model_from_java(env, model, jmodel);

	/* Mesh domain */
	vcn_mesh_t* mesh = vcn_mesh_create();
	vcn_mesh_set_geometric_constraint(mesh,
					  VCN_MESH_GEOM_CONSTRAINT_MAX_EDGE_LENGTH,
					  0.1);
	vcn_mesh_generate_from_model(mesh, model);
	vcn_msh3trg_t* msh3trg = 
		vcn_mesh_get_msh3trg(mesh, true, true, true, true, true);
	vcn_mesh_destroy(mesh);

	vcn_bcond_t* bmeshcond =
		vcn_fem_bcond_create_from_model_to_mesh(msh3trg, bconditions);

	/* FEM Analysis */
	vcn_fem_elem_t* elemtype = vcn_fem_elem_create(VCN_TRG_LINEAR);

	double* displacement = 
		malloc(msh3trg->N_vertices * 2 * sizeof(*displacement));
	double* strain = 
		malloc(msh3trg->N_triangles * 3 * sizeof(*strain));

	vcn_fem_compute_2D_Solid_Mechanics(msh3trg, elemtype, material,
					   bmeshcond, false, NULL,
					   vcn_sparse_solve_Cholesky, 1,
					   jthickness,
					   2, NULL,
					   displacement,
					   strain,
					   "static_elasticity2D_UT.log");

	jobject jmesh_results = jMeshResults_new(env);
	load_jMesh_from_msh3trg(env, msh3trg, jmesh_results);
	set_results_into_jmesh(env, jmesh_results, displacement, strain);

	/* Free memory */
	vcn_fem_bcond_destroy(bconditions);
	vcn_fem_bcond_destroy(bmeshcond);
	vcn_model_destroy(model);
	vcn_msh3trg_destroy(msh3trg);
	vcn_fem_material_destroy(material);
	vcn_fem_elem_destroy(elemtype);
	free(displacement);
	free(strain);
	return jmesh_results;
}

static void get_bconditions_from_java(JNIEnv *env, vcn_bcond_t *bcond,
				      jobject jbCondVtx, jobject jbCondSgm)
{
	/* Dirichlet on vertices */
	bcond->N_Dirichlet_on_vtx = jBoundaryConditions_getN(env, jbCondVtx);
	jintArray jbc_idx = jBoundaryConditions_getIdsRef(env, jbCondVtx);
	jint *bc_idx = (*env)->GetIntArrayElements(env, jbc_idx, NULL);
	jcharArray jbc_dof = jBoundaryConditions_getDofRef(env, jbCondVtx);
	jchar *bc_dof = (*env)->GetCharArrayElements(env, jbc_dof, NULL);
	jdoubleArray jbc_val = jBoundaryConditions_getValuesRef(env, jbCondVtx);
	jdouble *bc_val = (*env)->GetDoubleArrayElements(env, jbc_val, NULL);
	for (int i = 0; i < bcond->N_Dirichlet_on_vtx; i++) {
		bcond->Dirichlet_on_vtx_idx[i] = bc_idx[i];
		set_bc_dof(&(bcond->Dirichlet_on_vtx_dof_mask[i * 2]), bc_dof[i]);
		bcond->Dirichlet_on_vtx_val[i * 2] = bc_val[i * 2];
		bcond->Dirichlet_on_vtx_val[i*2+1] = bc_val[i*2+1];
	}
	(*env)->ReleaseIntArrayElements(env, jbc_idx, bc_idx, JNI_ABORT);
	(*env)->ReleaseCharArrayElements(env, jbc_dof, bc_dof, JNI_ABORT);
	(*env)->ReleaseDoubleArrayElements(env, jbc_val, bc_val, JNI_ABORT);

	/* Dirichlet on Segments */
	bcond->N_Dirichlet_on_sgm = jBoundaryConditions_getN(env, jbCondSgm);
	jbc_idx = jBoundaryConditions_getIdsRef(env, jbCondSgm);
	bc_idx = (*env)->GetIntArrayElements(env, jbc_idx, NULL);
	jbc_dof = jBoundaryConditions_getDofRef(env, jbCondSgm);
	jbc_val = jBoundaryConditions_getValuesRef(env, jbCondSgm);
	bc_val = (*env)->GetDoubleArrayElements(env, jbc_val, NULL);
	for (int i = 0; i < bcond->N_Dirichlet_on_sgm; i++) {
		bcond->Dirichlet_on_sgm_idx[i] = bc_idx[i];
		set_bc_dof(&(bcond->Dirichlet_on_sgm_dof_mask[i * 2]),
			   bc_dof[i]);
		bcond->Dirichlet_on_sgm_val[i * 2] = bc_val[i * 2];
		bcond->Dirichlet_on_sgm_val[i*2+1] = bc_val[i*2+1];
	}
	(*env)->ReleaseIntArrayElements(env, jbc_idx, bc_idx, JNI_ABORT);
	(*env)->ReleaseCharArrayElements(env, jbc_dof, bc_dof, JNI_ABORT);
	(*env)->ReleaseDoubleArrayElements(env, jbc_val, bc_val, JNI_ABORT);
}

static jint jBoundaryConditions_getN(JNIEnv *env, jobject jBoundaryConditions)
{
	jclass class = (*env)->GetObjectClass(env, jBoundaryConditions);
	jmethodID method_id = (*env)->GetMethodID(env, class, "getN", "()I");
	jint N = (*env)->CallIntMethod(env, jBoundaryConditions, method_id);
	return N;
}

static jintArray jBoundaryConditions_getIdsRef(JNIEnv *env,
					       jobject jBoundaryConditions)
{
	jclass class = (*env)->GetObjectClass(env, jBoundaryConditions);
	jmethodID method_id = 
		(*env)->GetMethodID(env, class, "getIdsRef", "()[I");
	jintArray idsRef =
		(*env)->CallObjectMethod(env, jBoundaryConditions, method_id);
	return idsRef;
}

static jcharArray jBoundaryConditions_getDofRef(JNIEnv *env,
						jobject jBoundaryConditions)
{
	jclass class = (*env)->GetObjectClass(env, jBoundaryConditions);
	jmethodID method_id =
		(*env)->GetMethodID(env, class, "getDofRef", "()[C");
	jcharArray dofRef =
		(*env)->CallObjectMethod(env, jBoundaryConditions, method_id);
	return dofRef;
}

static jdoubleArray jBoundaryConditions_getValuesRef(JNIEnv *env,
						     jobject jBoundaryConditions)
{
	jclass class = (*env)->GetObjectClass(env, jBoundaryConditions);
	jmethodID method_id =
		(*env)->GetMethodID(env, class, "getValuesRef", "()[D");
	jdoubleArray valuesRef =
		(*env)->CallObjectMethod(env, jBoundaryConditions, method_id);
	return valuesRef;
}

static void set_bc_dof(bool enabled[2], jchar dof)
{
	switch(dof) {
	case 'x':
		enabled[0] = true;
		enabled[1] = false;
		break;
	case 'y':
		enabled[0] = false;
		enabled[1] = true;
		break;
	case 'a':
		enabled[0] = true;
		enabled[1] = true;
		break;
	default:
		enabled[0] = false;
		enabled[1] = false;
	}	
}

static void get_material_from_java(JNIEnv *env, vcn_fem_material_t *material,
				   jobject jmaterial)
{
	vcn_fem_material_set_poisson_module
		(material, jMaterial_getPoissonModulus(env, jmaterial));
	vcn_fem_material_set_elasticity_module
		(material, jMaterial_getYoungModulus(env, jmaterial));
}

static jdouble jMaterial_getPoissonModulus(JNIEnv *env, jobject jmaterial)
{
	jclass class = (*env)->GetObjectClass(env, jmaterial);
	jmethodID method_id = 
		(*env)->GetMethodID(env, class, "getPoissonModulus", "()D");
	jdouble poisson_modulus =
		(*env)->CallDoubleMethod(env, jmaterial, method_id);
	return poisson_modulus;
}

static jdouble jMaterial_getYoungModulus(JNIEnv *env, jobject jmaterial)
{
	jclass class = (*env)->GetObjectClass(env, jmaterial);
	jmethodID method_id =
		(*env)->GetMethodID(env, class, "getYoungModulus", "()D");
	jdouble young_modulus =
		(*env)->CallDoubleMethod(env, jmaterial, method_id);
	return young_modulus;
}

static jobject jMeshResults_new(JNIEnv *env)
{
	jclass class =
		(*env)->FindClass(env, "Lvcn/pdeBot/finiteElement/MeshResults;");
	jmethodID method_id = (*env)->GetMethodID(env, class, "<init>", "()V");
	jobject instance = (*env)->NewObject(env, class, method_id);
	return instance;

}

static void set_results_into_jmesh(JNIEnv *env, jobject jmesh,
				   const double *displacement,
				   const double *strain)
{
	
}
