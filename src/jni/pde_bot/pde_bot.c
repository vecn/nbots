#include <stdlib.h>
#include "vcn/eigen_bot.h"
#include "vcn/geometric_bot.h"
#include "vcn/pde_bot.h"

#include "pde_bot.h"

static void get_bconditions_from_java(JNIEnv *env, vcn_bcond_t *bcond,
				      jobject jbCondVtx, jobject jbCondSgm);

static jint jBoundaryConditions_get_N(JNIEnv *env, jobject jBoundaryConditions);

static void get_material_from_java(JNIEnv *env, vcn_fem_material_t *material,
				   jobject jmaterial);

static void get_model_from_java(JNIEnv *env, vcn_model_t *model, jobject jmodel);

static jobject jMeshResults_new(JNIEnv *env);

static void cast_msh3trg_to_jmesh(JNIEnv *env, const vcn_msh3trg_t *const msh3trg,
				  jobject jmesh);
static void set_results_into_jmesh(JNIEnv *env, jobject jmesh, const double *displacement,
				   const double *strain);

/*
 * Class:     vcn_pdeBot_finiteElement_solidMechanics_StaticElasticity2D
 * Method:    solve
 * Signature: (Lvcn/geometricBot/Model;Lvcn/pdeBot/Material;Lvcn/pdeBot/BoundaryConditions;Lvcn/pdeBot/BoundaryConditions;D)Lvcn/pdeBot/finiteElement/MeshResults;
 */
JNIEXPORT jobject JNICALL Java_vcn_pdeBot_finiteElement_solidMechanics_StaticElasticity2D_solve
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
	cast_msh3trg_to_jmesh(env, msh3trg, jmesh_results);
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
	bcond->N_Dirichlet_on_vtx = jBoundaryConditions_get_N(env, jbCondVtx);
}

static jint jBoundaryConditions_get_N(JNIEnv *env, jobject jBoundaryConditions)
{
	jclass class = (*env)->GetObjectClass(env, jBoundaryConditions);
	jmethodID method_id = (*env)->GetMethodID(env, class, "getN", "()I");
	jint N = (*env)->CallIntMethod(env, jBoundaryConditions, method_id);
	return N;
}

static void get_material_from_java(JNIEnv *env, vcn_fem_material_t *material,
				   jobject jmaterial)
{

}

static void get_model_from_java(JNIEnv *env, vcn_model_t *model, jobject jmodel)
{

}

static jobject jMeshResults_new(JNIEnv *env)
{
	jclass class =
		(*env)->FindClass(env, "Lvcn/pdeBot/finiteElement/MeshResults;");
	jmethodID method_id = (*env)->GetMethodID(env, class, "<init>", "()V");
	jobject instance = (*env)->NewObject(env, class, method_id);
	return instance;

}

static void cast_msh3trg_to_jmesh(JNIEnv *env, const vcn_msh3trg_t *const msh3trg,
				  jobject jmesh)
{

}

static void set_results_into_jmesh(JNIEnv *env, jobject jmesh, const double *displacement,
				   const double *strain)
{

}
