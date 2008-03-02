#ifndef PROCRUSTES_H
#define PROCRUSTES_H
void procr_transform_wrapper(int *S, int *n, int *d, int *G, double *vZo,
			     double *vZ_mcmc, double *vZ_mu_mcmc, int *verbose);


void procr_alloc(int n, int d, int G,
		 double ***A, double ***tZ, double ***tZo, double ***Ahalf, 
		 double ***AhalfInv, double ***tptrans, double ***eAvectors, 
		 double ***eADvalues, double ***teAvectors, double **avZ, 
		 double **eAvalues, double ***dd_helper, 
		 double ***dn_helper, double ***dd2_helper,
		 double **workspace);

int procr_transform(double **Z, double **Z_mu, double **Zo, int n, int d, int G,
		    double **pZ, double **pZ_mu,
		    double **A, double **tZ, double **tZo, double **Ahalf, 
		    double **AhalfInv, double **tptrans, double **eAvectors, 
		    double **eADvalues, double **teAvectors, double *avZ, 
		    double *eAvalues, double **dd_helper, 
		    double **dn_helper, double **dd2_helper,
		    double *workspace);
#endif
