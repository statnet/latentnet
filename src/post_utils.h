#ifndef POST_UTILS_H
#define POST_UTILS_H
#include "ergmm_structs.h"
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
void post_pred_wrapper(int *S, 
		       int *n, int *p, int *d,
		       int *dir,
		       int *family, int *iconsts, double *dconsts,
		       double *vX,
		       double *Z_mcmc, 
		       double *coef_mcmc,
		       double *sender_mcmc, double *receiver_mcmc,
		       int *sociality,
		       int *vobserved_ties,
		       double *vEY,
		       int *verbose);
void klswitch_wrapper(int *maxit, int *S, int *n, int *d, int *G,
		      double *vZ_mcmc, int *Z_ref, double *vZ_mu_mcmc, double *vZ_var_mcmc,
		      int *vZ_K, double *vZ_pK,
		      double *vQ,
		      int *verbose);

void klswitch_step1(ERGMM_MCMC_Par *samples, int S, int n, int d, int G, double **Q, double ***pK);
int klswitch_step2(double **Q, ERGMM_MCMC_Par *samples, ERGMM_MCMC_Par *tmp, 
		   int S, int n, int d, int G,
		   unsigned int *perm, unsigned int *bestperm, unsigned int *dir,
		   double ***pK);
#endif
