void vcn_mesh_duplicate_one_point_connections(vcn_mesh_t* mesh)
{
	/* Allocate lists to store triangles per vertex */
	vcn_iter_t* iter = vcn_iter_create();
	vcn_iter_set_bins(iter, mesh->ug_vtx);
	while (vcn_iter_has_more(iter)) {
		msh_vtx_t* vtx = vcn_iter_get_next(iter);
		void** attr = malloc(2 * sizeof(*attr));
		vcn_container_t* trg_x_vtx = vcn_container_create(VCN_CONTAINER_QUEUE);
		attr[0] = trg_x_vtx;
		attr[1] = vtx->attr;
		vtx->attr = attr;
	}
	/* Iterate over triangles to found relations */
	vcn_iterator_t* ht_iter = vcn_iterator_create();
	vcn_iterator_set_container(ht_iter, mesh->ht_trg);
	while (vcn_iterator_has_more(ht_iter)) {
		msh_trg_t* trg = (msh_trg_t*)vcn_iterator_get_next(ht_iter);
		vcn_container_insert((vcn_container_t*)((void**)trg->v1->attr)[0], trg);
		vcn_container_insert((vcn_container_t*)((void**)trg->v2->attr)[0], trg);
		vcn_container_insert((vcn_container_t*)((void**)trg->v3->attr)[0], trg);
	}
	vcn_iterator_destroy(ht_iter);

	/* Detect one point connections and duplicate vertices */
	vcn_container_t* new_vertices = vcn_container_create(VCN_CONTAINER_QUEUE);
	vcn_iter_restart(iter);
	while (vcn_iter_has_more(iter)) {
		msh_vtx_t* vtx = vcn_iter_get_next(iter);
		vcn_container_t* trg_x_vtx = (vcn_container_t*)((void**)vtx->attr)[0];
		if(vcn_container_get_length(trg_x_vtx) < 2)
			continue;
		msh_trg_t* trg = vcn_container_get_first(trg_x_vtx);

		msh_trg_t* trg_twist = trg;
		bool twist_around = false;
		while (NULL != mtrg_get_right_triangle(trg_twist, vtx)) {
			trg_twist = mtrg_get_right_triangle(trg_twist, vtx);
			if (trg_twist == trg) {
				twist_around = true;
				break;
			}
		}
		if (trg_twist == trg && twist_around)
			continue;

		vcn_container_t* trg_fan = vcn_container_create(VCN_CONTAINER_QUEUE);
		do {
			vcn_container_insert(trg_fan, trg_twist);
			trg_twist = mtrg_get_left_triangle(trg_twist, vtx);
		} while(NULL != trg_twist);

		if (vcn_container_get_length(trg_fan) == vcn_container_get_length(trg_x_vtx)) {
			vcn_container_destroy(trg_fan);
			continue;
		}
    
		msh_vtx_t* new_vtx = vcn_point2D_create();
		memcpy(new_vtx->x, vtx->x, 2 * sizeof(*(vtx->x)));
		new_vtx->attr = ((void**)vtx->attr)[1];
		vcn_container_insert(new_vertices, new_vtx);
		
		while (vcn_container_is_not_empty(trg_fan)) {
			trg = vcn_container_delete_first(trg_fan);
    
			msh_edge_t* s1 = mtrg_get_left_edge(trg, vtx);
			msh_edge_t* s2 = mtrg_get_right_edge(trg, vtx);
    
			if (trg->v1 == vtx)
				trg->v1 = new_vtx;
			else if (trg->v2 == vtx)
				trg->v2 = new_vtx;
			else if (trg->v3 == vtx)
				trg->v3 = new_vtx;
    
			if (s1->v1 == vtx)
				s1->v1 = new_vtx;
			else if (s1->v2 == vtx)
				s1->v2 = new_vtx;
    
			if (s2->v1 == vtx)
				s2->v1 = new_vtx;
			else if (s2->v2 == vtx)
				s2->v2 = new_vtx;
		}
    		vcn_container_destroy(trg_fan);
	}

	/* Free memory */
	vcn_iter_restart(iter);
	while (vcn_iter_has_more(iter)) {
		msh_vtx_t* vtx = vcn_iter_get_next(iter);
		void** attr = vtx->attr;
		vcn_container_destroy(attr[0]);
		vtx->attr = attr[1];
		free(attr);
	}
	vcn_iter_destroy(iter);

	while (vcn_container_is_not_empty(new_vertices)) {
		msh_vtx_t* new_vtx = vcn_container_delete_first(new_vertices);    
		vcn_bins2D_insert(mesh->ug_vtx, new_vtx);
	}

	vcn_container_destroy(new_vertices);
}
