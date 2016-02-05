static double N1(double psi, double eta);
static double N2(double psi, double eta);
static double N3(double psi, double eta);
static double dN1_dpsi(double psi, double eta);
static double dN1_deta(double psi, double eta);
static double dN2_dpsi(double psi, double eta);
static double dN2_deta(double psi, double eta);
static double dN3_dpsi(double psi, double eta);
static double dN3_deta(double psi, double eta);

vcn_fem_elem_t* vcn_fem_elem_create_triangle(){
	/* Triangular element of linear interpolation */
	vcn_fem_elem_t* elemtype = (vcn_fem_elem_t*)
		malloc(sizeof(vcn_fem_elem_t));
	elemtype->N_nodes = 3;
	elemtype->N_Gauss_points = 1;

	elemtype->Ni = 
		(double (**)(double, double))calloc(elemtype->N_nodes, sizeof(void*));
	elemtype->dNi_dpsi = 
		(double (**)(double, double))calloc(elemtype->N_nodes, sizeof(void*));
	elemtype->dNi_deta = 
		(double (**)(double, double))calloc(elemtype->N_nodes, sizeof(void*));  
	elemtype->Ni[0] = N1;
	elemtype->Ni[1] = N2;
	elemtype->Ni[2] = N3;
	elemtype->dNi_dpsi[0] = dN1_dpsi;
	elemtype->dNi_dpsi[1] = dN2_dpsi;
	elemtype->dNi_dpsi[2] = dN3_dpsi;
	elemtype->dNi_deta[0] = dN1_deta;
	elemtype->dNi_deta[1] = dN2_deta;
	elemtype->dNi_deta[2] = dN3_deta;
  
	elemtype->psi = 
		(double*)malloc(elemtype->N_Gauss_points*sizeof(double));
	elemtype->eta = 
		(double*)malloc(elemtype->N_Gauss_points*sizeof(double));
	elemtype->gp_weight = 
		(double*)malloc(elemtype->N_Gauss_points*sizeof(double));

	elemtype->psi[0] = 1.0/3.0;
	elemtype->eta[0] = 1.0/3.0;

	elemtype->gp_weight[0] = 0.5;

	return elemtype;
}

static inline double N1(double psi, double eta){
	return 1.0 - psi - eta;
}

static inline double N2(double psi, double eta){
	return psi;
}

static inline double N3(double psi, double eta){
	return eta;
}

static inline double dN1_dpsi(double psi, double eta){
	return -1.0;
}

static inline double dN1_deta(double psi, double eta){
	return -1.0;
}

static inline double dN2_dpsi(double psi, double eta){
	return 1.0;
}

static inline double dN2_deta(double psi, double eta){
	return 0.0;
}

static inline double dN3_dpsi(double psi, double eta){
	return 0.0;
}

static inline double dN3_deta(double psi, double eta){
	return 1.0;
}
