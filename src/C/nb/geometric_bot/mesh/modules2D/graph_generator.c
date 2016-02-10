vcn_graph_t* vcn_mesh_create_vtx_graph(const vcn_mesh_t *const restrict mesh)
{
	vcn_graph_t* graph = (vcn_graph_t*) calloc(1, sizeof(vcn_graph_t));

	graph->N = vcn_bins2D_length(mesh->ug_vtx);
	graph->N_adj = calloc(graph->N, sizeof(*(graph->N_adj)));
	graph->adj = malloc(graph->N * sizeof(*(graph->adj)));
	vcn_iterator_t* iter = vcn_iterator_create();
	vcn_iterator_set_container(iter, mesh->ht_edge);
	while (vcn_iterator_has_more(iter)) {
		msh_edge_t* edge = (msh_edge_t*)vcn_iterator_get_next(iter);
		uint32_t idx1 = ((uint32_t*)((void**)edge->v1->attr)[0])[0];
		uint32_t idx2 = ((uint32_t*)((void**)edge->v2->attr)[0])[0];
		graph->N_adj[idx1] += 1;
		graph->N_adj[idx2] += 1;
	}

	for (uint32_t i = 0; i < graph->N; i++)
		graph->adj[i] = malloc(graph->N_adj[i] * sizeof(*(graph->adj[i])));

	uint32_t* adj_next_idx = calloc(graph->N, sizeof(*adj_next_idx));
	vcn_iterator_restart(iter);
	while (vcn_iterator_has_more(iter)) {
		msh_edge_t* edge = (msh_edge_t*)vcn_iterator_get_next(iter);
		uint32_t idx1 = ((uint32_t*)((void**)edge->v1->attr)[0])[0];
		uint32_t idx2 = ((uint32_t*)((void**)edge->v2->attr)[0])[0];
		graph->adj[idx1][adj_next_idx[idx1]] = idx2;
		graph->adj[idx2][adj_next_idx[idx2]] = idx1;
		adj_next_idx[idx1] += 1;
		adj_next_idx[idx2] += 1;
	}
	vcn_iterator_destroy(iter);
	free(adj_next_idx);

	return graph;
}

vcn_graph_t* vcn_mesh_create_elem_graph(const vcn_mesh_t *const restrict mesh)
{
	vcn_graph_t* graph = calloc(1, sizeof(*graph));
	graph->N = vcn_container_get_length(mesh->ht_trg);
	graph->N_adj = calloc(graph->N, sizeof(*(graph->N_adj)));
	graph->adj = malloc(graph->N * sizeof(*(graph->adj)));
	vcn_iterator_t* iter = vcn_iterator_create();
	vcn_iterator_set_container(iter, mesh->ht_trg);
	while (vcn_iterator_has_more(iter)) {
		msh_trg_t* trg = (msh_trg_t*)vcn_iterator_get_next(iter);
		uint32_t id = ((uint32_t*)((void**)trg->attr)[0])[0];
		if (NULL != trg->t1)
			graph->N_adj[id] += 1;
		if (NULL != trg->t2)
			graph->N_adj[id] += 1;
		if (NULL != trg->t3)
			graph->N_adj[id] += 1;
	}

	vcn_iterator_restart(iter);
	while (vcn_iterator_has_more(iter)) {
		msh_trg_t* trg = (msh_trg_t*)vcn_iterator_get_next(iter);
		uint32_t id = ((uint32_t*)((void**)trg->attr)[0])[0];
		graph->adj[id] = malloc(graph->N_adj[id] * sizeof(*(graph->adj[id])));
		int cnt = 0;
		if (NULL != trg->t1) {
			uint32_t id2 = ((uint32_t*)((void**)trg->t1->attr)[0])[0];
			graph->adj[id][cnt++] = id2;
		}
		if (NULL != trg->t2) {
			uint32_t id2 = ((uint32_t*)((void**)trg->t2->attr)[0])[0];
			graph->adj[id][cnt++] = id2;
		}
		if (NULL != trg->t3){
			uint32_t id2 = ((uint32_t*)((void**)trg->t3->attr)[0])[0];
			graph->adj[id][cnt++] = id2;
		}
	}
	vcn_iterator_destroy(iter);
	return graph;
}
