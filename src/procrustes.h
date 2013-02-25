#ifndef PROCRUSTES_H
#define PROCRUSTES_H
void procr_transform_wrapper(int *S, int *n, int *d, int *G, double *vZo,
			     double *vZ_mcmc, double *vZ_mu_mcmc, int *verbose);


void procr_alloc(int n, int d, int G,
		 double **avZ, double ***tZZo, double ***U, double **S, double ***tV,
		 double **workspace);

int procr_transform(double **Z, double **Z_mean, double **Zo, int n, int d, int G,
		    double **pZ, double **pZ_mean,
		    double *avZ, double **tZZo, double **U, double *S, double **tV,
		    double *workspace);
#endif
