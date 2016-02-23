#include <stdlib.h>
#include <alloca.h>

#include "nb/geometric_bot.h"
#include "nb/pde_bot.h"

#include "geometric_bot/jMesh.h"
#include "geometric_bot/jModel.h"
#include "pde_bot/JNI_jPdeBot.h"

static void get_bconditions_from_java(JNIEnv *env, nb_bcond_t *bcond,
				      jobject jBCDirichletVtx,
				      jobject jBCNeumannVtx,
				      jobject jBCDirichletSgm,
				      jobject jBCNeumannSgm);
static void set_jBoundaryCondition_into_bcond(JNIEnv *env, nb_bcond_t *bcond,
					      jobject jBCDirichletVtx,
					      nb_bcond_id bc_type,
					      nb_bcond_where bc_where);

static jint jBoundaryConditions_getN(JNIEnv *env, jobject jBoundaryConditions);
static jintArray jBoundaryConditions_getIdsRef(JNIEnv *env,
					       jobject jBoundaryConditions);
static jcharArray jBoundaryConditions_getDofRef
			(JNIEnv *env, jobject jBoundaryConditions);
static jdoubleArray jBoundaryConditions_getValuesRef
			(JNIEnv *env, jobject jBoundaryConditions);
static void set_bc_dof(bool enabled[2], jchar dof);
static void get_material_from_java(JNIEnv *env, vcn_fem_material_t *material,
				   jobject jMaterial);
static jdouble jMaterial_getPoissonModulus(JNIEnv *env, jobject jMaterial);
static jdouble jMaterial_getYoungModulus(JNIEnv *env, jobject jMaterial);
static jobject jMeshResults_create(JNIEnv *env);
static void set_results_into_jMesh(JNIEnv *env, jobject jMesh,
				   const double *displacement,
				   const double *strain, jint N_vtx, jint N_trg);
static void load_displacements_into_jMeshResults(JNIEnv *env, jobject jMesh,
						 const double *displacement,
						 jint N);
static void load_strain_into_jMeshResults(JNIEnv *env, jobject jMesh,
					  const double *strain,
					  jint N);

/*
 * Class:     vcn_pdeBot_finiteElement_solidMechanics_StaticElasticity2D
 * Method:    solve
 * Signature: (Lnb/geometricBot/Model;Lnb/pdeBot/Material;Lnb/pdeBot/BoundaryConditions;Lnb/pdeBot/BoundaryConditions;Lnb/pdeBot/BoundaryConditions;Lnb/pdeBot/BoundaryConditions;DI)Lnb/pdeBot/finiteElement/MeshResults;
 */
JNIEXPORT jobject JNICALL
Java_nb_pdeBot_finiteElement_solidMechanics_StaticElasticity2D_solve
		(JNIEnv *env, jclass class, jobject jModel, jobject jMaterial,
		 jobject jBCDirichletVtx, jobject jBCNeumannVtx,
		 jobject jBCDirichletSgm, jobject jBCNeumannSgm,
		 jdouble jthickness, jint jN_nodes)
{
	/* Read optimization parameters */
	uint16_t bcond_size = nb_bcond_get_memsize(2);
	nb_bcond_t* bcond = alloca(bcond_size);
	nb_bcond_init(bcond, 2);
	get_bconditions_from_java(env, bcond,
				  jBCDirichletVtx, jBCNeumannVtx,
				  jBCDirichletSgm, jBCNeumannSgm);
	
	vcn_fem_material_t* material = vcn_fem_material_create();
	get_material_from_java(env, material, jMaterial);
	
	vcn_model_t* model = vcn_model_create();
	vcn_model_load_from_jModel(env, model, jModel);

	/* Mesh domain */
	vcn_mesh_t* mesh = vcn_mesh_create();
	if (0 < jN_nodes) {
		vcn_mesh_set_size_constraint(mesh,
					     NB_MESH_SIZE_CONSTRAINT_MAX_VTX,
					     jN_nodes);
		vcn_mesh_set_geometric_constraint(mesh,
						  NB_MESH_GEOM_CONSTRAINT_MAX_EDGE_LENGTH,
						  NB_GEOMETRIC_TOL);
	}
	vcn_mesh_generate_from_model(mesh, model);
	vcn_msh3trg_t* msh3trg = 
		vcn_mesh_get_msh3trg(mesh, true, true, true, true, true);
	vcn_mesh_destroy(mesh);

	/* FEM Analysis */
	vcn_fem_elem_t* elemtype = vcn_fem_elem_create(NB_TRG_LINEAR);

	double* displacement = 
		malloc(msh3trg->N_vertices * 2 * sizeof(*displacement));
	double* strain = 
		malloc(msh3trg->N_triangles * 3 * sizeof(*strain));

	vcn_fem_compute_2D_Solid_Mechanics(msh3trg, elemtype, material,
					   bcond, false, NULL, true,
					   jthickness, NULL,
					   displacement, strain);

	jobject jMeshResults = jMeshResults_create(env);
	load_jMesh_from_msh3trg(env, msh3trg, jMeshResults);
	set_results_into_jMesh(env, jMeshResults, displacement, strain,
			       msh3trg->N_vertices, msh3trg->N_triangles);

	nb_bcond_finish(bcond);

	vcn_model_destroy(model);
	vcn_msh3trg_destroy(msh3trg);
	vcn_fem_material_destroy(material);
	vcn_fem_elem_destroy(elemtype);
	free(displacement);
	free(strain);
	return jMeshResults;
}

static void get_bconditions_from_java(JNIEnv *env, nb_bcond_t *bcond,
				      jobject jBCDirichletVtx,
				      jobject jBCNeumannVtx,
				      jobject jBCDirichletSgm,
				      jobject jBCNeumannSgm)
{
	set_jBoundaryCondition_into_bcond(env, bcond, jBCDirichletVtx,
					  NB_DIRICHLET, NB_BC_ON_POINT);
	set_jBoundaryCondition_into_bcond(env, bcond, jBCNeumannVtx,
					  NB_NEUMANN, NB_BC_ON_POINT);
	set_jBoundaryCondition_into_bcond(env, bcond, jBCDirichletSgm,
					  NB_DIRICHLET, NB_BC_ON_SEGMENT);
	set_jBoundaryCondition_into_bcond(env, bcond, jBCNeumannSgm,
					  NB_NEUMANN, NB_BC_ON_SEGMENT);
}

static void set_jBoundaryCondition_into_bcond(JNIEnv *env, nb_bcond_t *bcond,
					      jobject jBoundaryCondition,
					      nb_bcond_id bc_type,
					      nb_bcond_where bc_where)
{
	uint32_t N = jBoundaryConditions_getN(env, jBoundaryCondition);
	if (0 < N) {
		jintArray jbc_idx =
			jBoundaryConditions_getIdsRef(env, jBoundaryCondition);
		jint *bc_idx =
			(*env)->GetIntArrayElements(env, jbc_idx, NULL);
		jcharArray jbc_dof = 
			jBoundaryConditions_getDofRef(env, jBoundaryCondition);
		jchar *bc_dof =
			(*env)->GetCharArrayElements(env, jbc_dof, NULL);
		jdoubleArray jbc_val =
			jBoundaryConditions_getValuesRef(env,
							 jBoundaryCondition);
		jdouble *bc_val =
			(*env)->GetDoubleArrayElements(env, jbc_val, NULL);
		
		for (int i = 0; i < N; i++) {
			bool mask[2];
			set_bc_dof(mask, bc_dof[i]);
			double val[2];
			val[0] = bc_val[i * 2];
			val[1] = bc_val[i*2+1];
			nb_bcond_push(bcond, bc_type, bc_where,
				      bc_idx[i], mask, val);
		}
		(*env)->ReleaseIntArrayElements(env, jbc_idx,
						bc_idx, JNI_ABORT);
		(*env)->ReleaseCharArrayElements(env, jbc_dof,
						 bc_dof, JNI_ABORT);
		(*env)->ReleaseDoubleArrayElements(env, jbc_val,
						   bc_val, JNI_ABORT);
	}
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
				   jobject jMaterial)
{
	vcn_fem_material_set_poisson_module
		(material, jMaterial_getPoissonModulus(env, jMaterial));
	vcn_fem_material_set_elasticity_module
		(material, jMaterial_getYoungModulus(env, jMaterial));
}

static jdouble jMaterial_getPoissonModulus(JNIEnv *env, jobject jMaterial)
{
	jclass class = (*env)->GetObjectClass(env, jMaterial);
	jmethodID method_id = 
		(*env)->GetMethodID(env, class, "getPoissonModulus", "()D");
	jdouble poisson_modulus =
		(*env)->CallDoubleMethod(env, jMaterial, method_id);
	return poisson_modulus;
}

static jdouble jMaterial_getYoungModulus(JNIEnv *env, jobject jMaterial)
{
	jclass class = (*env)->GetObjectClass(env, jMaterial);
	jmethodID method_id =
		(*env)->GetMethodID(env, class, "getYoungModulus", "()D");
	jdouble young_modulus =
		(*env)->CallDoubleMethod(env, jMaterial, method_id);
	return young_modulus;
}

static jobject jMeshResults_create(JNIEnv *env)
{
	jclass class =
		(*env)->FindClass(env, "nb/pdeBot/finiteElement/MeshResults");
	jmethodID method_id = (*env)->GetMethodID(env, class, "<init>", "()V");
	jobject instance = (*env)->NewObject(env, class, method_id);
	return instance;

}

static void set_results_into_jMesh(JNIEnv *env, jobject jMesh,
				   const double *displacement,
				   const double *strain, jint N_vtx, jint N_trg)
{
	load_displacements_into_jMeshResults(env, jMesh, displacement, N_vtx);
	load_strain_into_jMeshResults(env, jMesh, strain, N_trg);
}

static void load_displacements_into_jMeshResults(JNIEnv *env, jobject jMesh,
						 const double *displacement,
						 jint N)
{
	jclass class =
		(*env)->GetObjectClass(env, jMesh);
	jfieldID field_id = 
		(*env)->GetFieldID(env, class, "displacement", "[F");
	jfloatArray jdisplacement =
		(*env)->NewFloatArray(env, 2 * N);
	jfloat *fdisp =
		(*env)->GetFloatArrayElements(env, jdisplacement, NULL);

	for (uint32_t i = 0; i < 2 * N; i++)
		fdisp[i] = displacement[i];
	(*env)->ReleaseFloatArrayElements(env, jdisplacement, fdisp, 0);
	(*env)->SetObjectField(env, jMesh, field_id, jdisplacement);
}

static void load_strain_into_jMeshResults(JNIEnv *env, jobject jMesh,
					  const double *strain,
					  jint N)
{
	jclass class =
		(*env)->GetObjectClass(env, jMesh);
	jfieldID field_id = 
		(*env)->GetFieldID(env, class, "strain", "[F");
	jfloatArray jstrain =
		(*env)->NewFloatArray(env, 3 * N);
	jfloat *fstrain =
		(*env)->GetFloatArrayElements(env, jstrain, NULL);

	for (uint32_t i = 0; i < 3 * N; i++)
		fstrain[i] = strain[i];
	(*env)->ReleaseFloatArrayElements(env, jstrain, fstrain, 0);
	(*env)->SetObjectField(env, jMesh, field_id, jstrain);
}
