#include <stdlib.h>
#include <string.h>
#include <alloca.h>

#include "nb/geometric_bot.h"
#include "nb/pde_bot.h"

#include "geometric_bot/jMesh.h"
#include "geometric_bot/jModel.h"
#include "pde_bot/JNI_jPdeBot.h"

typedef struct {
	double *disp;
	double *strain;
	double *stress;
	double *strain_on_nodes,
		*stress_on_nodes; /* Same memory buffer */
	double *von_mises;
} results_t;

static void results_init(results_t *results,
			 const vcn_msh3trg_t *msh3trg);
static void results_finish(results_t *results);

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
static void nb_analysis2D_load_from_jAnalysis2D(JNIEnv *env, 
						nb_analysis2D_t * analysis2D,
						nb_analysis2D_params *params2D,
						jobject jAnalysis2D);
static jdouble jMaterial_getPoissonModulus(JNIEnv *env, jobject jMaterial);
static jdouble jMaterial_getYoungModulus(JNIEnv *env, jobject jMaterial);
static jobject jMeshResults_create(JNIEnv *env);
static void jMeshResults_set_status(JNIEnv *env, jobject jMeshResults,
				    int fem_status);
static void set_results_into_jMesh(JNIEnv *env, jobject jMesh,
				   const double *displacement,
				   const double *strain,
				   const double *von_mises,
				   jint N_vtx);
static void load_displacements_into_jMeshResults(JNIEnv *env, jobject jMesh,
						 const double *displacement,
						 jint N);
static void load_strain_into_jMeshResults(JNIEnv *env, jobject jMesh,
					  const double *strain,
					  jint N);
static void load_von_mises_into_jMeshResults(JNIEnv *env, jobject jMesh,
					     const double *von_mises,
					     jint N);

/*
 * Class:     nb_pdeBot_finiteElement_solidMechanics_StaticElasticity2D
 * Method:    solve
 * Signature: (Lnb/geometricBot/Model;Lnb/pdeBot/Material;Lnb/pdeBot/finiteElement/solidMechanics/Analysis2D;Lnb/pdeBot/BoundaryConditions;Lnb/pdeBot/BoundaryConditions;Lnb/pdeBot/BoundaryConditions;Lnb/pdeBot/BoundaryConditions;DI)Lnb/pdeBot/finiteElement/MeshResults;
 */
JNIEXPORT jobject JNICALL
Java_nb_pdeBot_finiteElement_solidMechanics_StaticElasticity2D_solve
		(JNIEnv *env, jclass class, jobject jModel,
		 jobject jMaterial, jobject jAnalysis2D,
		 jobject jBCDirichletVtx, jobject jBCNeumannVtx,
		 jobject jBCDirichletSgm, jobject jBCNeumannSgm,
		 jint jN_nodes)
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
	
	vcn_model_t* model = alloca(vcn_model_get_memsize());
	vcn_model_init(model);
	vcn_model_load_from_jModel(env, model, jModel);

	nb_analysis2D_t analysis2D;
	nb_analysis2D_params params2D;
	nb_analysis2D_load_from_jAnalysis2D(env, &analysis2D, &params2D, jAnalysis2D);

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

	/* Allocate memory */
	results_t results;
	results_init(&results, msh3trg);

	int status_fem =
	  vcn_fem_compute_2D_Solid_Mechanics(msh3trg, elemtype, material,
					     bcond, false, NULL, analysis2D,
					     &params2D, NULL,
					     results.disp, results.strain);
	

	jobject jMeshResults = jMeshResults_create(env);
	jMeshResults_set_status(env, jMeshResults, status_fem);
	load_jMesh_from_msh3trg(env, msh3trg, jMeshResults);

	if (0 == status_fem) {
		vcn_fem_compute_stress_from_strain(msh3trg->N_triangles,
						   msh3trg->vertices_forming_triangles,
						   elemtype, material,
						   analysis2D, results.strain, NULL, 
						   results.stress);

		vcn_fem_interpolate_from_Gauss_points_to_nodes(msh3trg, elemtype,
							       3, results.stress,
							       results.stress_on_nodes);

		vcn_fem_compute_von_mises(msh3trg->N_vertices, results.stress_on_nodes,
					  results.von_mises);
	
		vcn_fem_interpolate_from_Gauss_points_to_nodes(msh3trg, elemtype,
							       3, results.strain,
							       results.strain_on_nodes);

		set_results_into_jMesh(env, jMeshResults,
				       results.disp, results.strain_on_nodes,
				       results.von_mises, msh3trg->N_vertices);
	}

	nb_bcond_finish(bcond);
	vcn_model_finish(model);
	vcn_msh3trg_destroy(msh3trg);
	vcn_fem_material_destroy(material);
	vcn_fem_elem_destroy(elemtype);
	results_finish(&results);
	return jMeshResults;
}

static void results_init(results_t *results,
			 const vcn_msh3trg_t *msh3trg)
{
	uint32_t size_disp = msh3trg->N_vertices * 2 *
		sizeof(*(results->disp));
	uint32_t size_strain = msh3trg->N_triangles * 3 *
		sizeof(*(results->strain));
	uint32_t size_stress = msh3trg->N_triangles * 3 *
		sizeof(*(results->stress));
	uint32_t size_data_on_nodes = msh3trg->N_vertices * 3 *
		sizeof(*(results->strain_on_nodes));
	uint32_t size_von_mises = msh3trg->N_vertices *
		sizeof(*(results->von_mises));
	uint32_t total_size = size_disp + size_strain +
		size_stress + size_data_on_nodes + size_von_mises;
	char *memblock = malloc(total_size);

	results->disp = (void*) memblock;
	results->strain = (void*)(memblock + size_disp);
	results->stress = (void*)(memblock + size_disp + size_strain);
	results->strain_on_nodes = (void*)(memblock + size_disp +
					   size_strain + size_stress);
	results->stress_on_nodes = (void*)(memblock + size_disp +
					   size_strain + size_stress);
	results->von_mises = (void*)(memblock + size_disp + size_strain +
				     size_stress + size_data_on_nodes);
}

static inline void results_finish(results_t *results)
{
	free(results->disp);
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

static void nb_analysis2D_load_from_jAnalysis2D(JNIEnv *env, 
						nb_analysis2D_t * analysis2D,
						nb_analysis2D_params *params2D,
						jobject jAnalysis2D)
{
	const char *str_class =
		"nb/pdeBot/finiteElement/solidMechanics/Analysis2D";
	jclass enum_class = (*env)->FindClass(env, str_class);
	jmethodID enum_get_name_id = 
		(*env)->GetMethodID(env, enum_class, "name",
				    "()Ljava/lang/String;");
	jstring enum_type = (*env)->CallObjectMethod(env, jAnalysis2D,
						     enum_get_name_id);
	const char *str_enum = (*env)->GetStringUTFChars(env, enum_type, 0);
	
	if (0 == strcmp("PLANE_STRESS", str_enum)) {
		*analysis2D = NB_PLANE_STRESS;
		jmethodID methodID = 
			(*env)->GetMethodID(env, enum_class,
					    "getThickness", "()D");
		jdouble jthickness =
			(*env)->CallDoubleMethod(env, jAnalysis2D, methodID);
		params2D->thickness = jthickness;
	} else if (0 == strcmp("PLANE_STRAIN", str_enum)) {
		*analysis2D = NB_PLANE_STRAIN;
		params2D->thickness = 1.0;
	} else if (0 == strcmp("SOLID_OF_REVOLUTION", str_enum)) {
		*analysis2D = NB_SOLID_OF_REVOLUTION;
		params2D->thickness = 1.0;
	} else {
		*analysis2D = NB_PLANE_STRESS;
		params2D->thickness = 1.0;
	}
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

static void jMeshResults_set_status(JNIEnv *env, jobject jMeshResults,
				    int fem_status)
{
	jclass class = (*env)->GetObjectClass(env, jMeshResults);
	jmethodID method_id =
		(*env)->GetMethodID(env, class, "setStatus", "(I)V");
	(*env)->CallVoidMethod(env, jMeshResults, method_id, fem_status);
}

static void set_results_into_jMesh(JNIEnv *env, jobject jMesh,
				   const double *displacement,
				   const double *strain,
				   const double *von_mises,
				   jint N_vtx)
{
	load_displacements_into_jMeshResults(env, jMesh, displacement, N_vtx);
	load_strain_into_jMeshResults(env, jMesh, strain, N_vtx);
	load_von_mises_into_jMeshResults(env, jMesh, von_mises, N_vtx);
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

static void load_von_mises_into_jMeshResults(JNIEnv *env, jobject jMesh,
					     const double *von_mises,
					     jint N)
{
	jclass class =
		(*env)->GetObjectClass(env, jMesh);
	jfieldID field_id = 
		(*env)->GetFieldID(env, class, "VonMisesStress", "[F");
	jfloatArray jstress =
		(*env)->NewFloatArray(env, N);
	jfloat *fstress =
		(*env)->GetFloatArrayElements(env, jstress, NULL);

	for (uint32_t i = 0; i < N; i++)
		fstress[i] = von_mises[i];
	(*env)->ReleaseFloatArrayElements(env, jstress, fstress, 0);
	(*env)->SetObjectField(env, jMesh, field_id, jstress);
}
