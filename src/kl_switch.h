/*  File src/kl_switch.h in package latentnet, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2017 Statnet Commons
 */
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
