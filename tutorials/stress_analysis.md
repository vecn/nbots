# Tutorial: Customized stress analysis with FEM

This tutorial explains the basics for running an stress analysis with FEM.

Use the function
~~~C
nb_mesh2D_t *mesh2D;
nb_fem_elem_t* elem = nb_fem_elem_create(NB_TRG_LINEAR);
nb_bcond_t *bcond;
nb_material_t* material = nb_material_create();
nb_analysis2D_t analysis2D;
nb_analysis2D_params params2D;

double *disp;
double *strain;
nb_fem_compute_2D_Solid_Mechanics(mesh2D, elem,
				  material, bcond,
				  false, NULL,
				  analysis2D,
				  &params2D, NULL,
				  disp, strain);
~~~

(Pending to develop)