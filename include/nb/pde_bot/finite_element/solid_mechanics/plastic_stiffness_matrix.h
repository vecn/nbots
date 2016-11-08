#ifndef PLASTIC_STIFFNESS_MATRIX_H_INCLUDED
#define PLASTIC_STIFFNESS_MATRIX_H_INCLUDED

int updated_stiffness_matrix(nb_sparse_t* K, nb_plastified_analysis2D *elem_regime,
                        const bool *elements_enabled, const nb_mesh2D_t *part,
                        const nb_fem_elem_t *elem,
                        const nb_material_t *material,
                        nb_analysis2D_t analysis2D,
                        nb_analysis2D_params *params2D,
                        double *elastic_strain,
                        uint32_t *simultaneous_elements,
                        double *aux_strain,
                        uint32_t *N_simultaneous_plastic_elem,
                        uint32_t *plastified_elem,
                        uint32_t *N_plastic_elem);
int updated_plastified_stiffness_matrix(nb_sparse_t* K, nb_plastified_analysis2D el_regime,
                                    bool is_enabled, const nb_mesh2D_t *part,
                                    uint32_t plastified_elem,
                                    const nb_fem_elem_t *elem,
                                    const nb_material_t *material,
                                    nb_analysis2D_t analysis2D,
                                    nb_analysis2D_params *params2D);
static int assemble_plastified_element(const nb_fem_elem_t *elem, uint32_t id,
                                        const nb_mesh2D_t *part,
                                        const nb_material_t *material,
                                        bool is_enabled,
                                        nb_plastified_analysis2D el_regime,
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
