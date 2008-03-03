#ifndef KL_SWITCH_H
#define KL_SWITCH_H
#include "ergmm_structs.h"

void klswitch_wrapper(int *maxit, int *S, int *n, int *d, int *G,
		      double *vZ_mcmc, int *Z_ref, double *vZ_mu_mcmc, double *vZ_var_mcmc,
		      int *vZ_K, double *vZ_pK,
		      double *vQ,
		      int *verbose);

void klswitch_step1(ERGMM_MCMC_Par *sample, int S, int n, int G, double **Q, double ***pK);
int klswitch_step2(double **Q, ERGMM_MCMC_Par *sample, ERGMM_MCMC_Par *tmp, 
		   unsigned int S, unsigned int n, unsigned int d, unsigned int G,
		   unsigned int *perm, unsigned int *bestperm, unsigned int *dir,
		   double ***pK);
#endif
