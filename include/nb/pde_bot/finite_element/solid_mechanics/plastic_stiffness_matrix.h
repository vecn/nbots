#ifndef PLASTIC_STIFFNESS_MATRIX_H_INCLUDED
#define PLASTIC_STIFFNESS_MATRIX_H_INCLUDED


static int assemble_plastified_element(const nb_fem_elem_t *elem, uint32_t id,
			    const nb_mesh2D_t *part,
			    const nb_material_t *material,
			    bool is_enabled,
			    nb_plastified_analysis2D elem_regime,
			    nb_analysis2D_t analysis2D,
			    nb_analysis2D_params *params2D,
			    nb_sparse_t *K);
static int integrate_plastic_elemental_system
            (const nb_fem_elem_t *elem, uint32_t id,
			 double D[4], const nb_mesh2D_t *part,
			 nb_analysis2D_params *params2D,
			 double *Ke);
void modify_global_system(const nb_fem_elem_t *elem, uint32_t id, const nb_mesh2D_t *part,
                          double *Ke_elastic, double *Ke_plastic, nb_sparse_t *K);


#endif // PLASTIC_STIFFNESS_MATRIX_H_INCLUDED
