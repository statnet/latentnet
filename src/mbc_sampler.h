/*  File src/mbc_sampler.h in package latentnet, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2014 Statnet Commons
 */
#ifndef MBC_SAMPLER_H
#define MBC_SAMPLER_H

#include "ergmm_structs.h"

/* First 1 positions in the outlists are reserved for special values:
   [0] Iteration with the highest likelihood so far.
*/
#define MBC_OUTLISTS_RESERVE 1
void MBC_MCMC_wrapper(int *sample_size,
		      int *interval,

		      int *n,
		      int *d,
		      int *G,
			  
		      double *lpZ_mcmc, 
		      double *lpLV_mcmc, 
		      
		      double *vZ,
		      
		      double *Z_pK,
		      double *Z_mean,
		      double *Z_var,
		      int *Z_K,
		      
		      double *Z_var_prior, 
		      double *Z_mean_prior_var, 
		      double *Z_K_prior,
		      double *Z_var_df,
		      
		      int *Z_K_mcmc, 
		      double *Z_pK_mcmc, 
		      double *Z_mean_mcmc, 
		      double *Z_var_mcmc);
void MBC_MCMC_init(unsigned int sample_size, 
		   unsigned int interval, 

		   unsigned int n,
		   unsigned int d,
		   unsigned int G,
		   
		   double *lpZ_mcmc,
		   double *lpLV_mcmc, 
		   
		   double **Z,

		   double *Z_pK, 
		   double **Z_mean_start, 
		   double *Z_var, 
		   unsigned int *Z_K,

		   double Z_var_prior,
		   double Z_mean_prior_var,
		   double Z_K_prior,
		   double Z_var_df,

		   int *K_mcmc,
		   double *Z_pK_mcmc,
		   double *Z_mean_mcmc,
		   double *Z_var_mcmc);

void MBC_MCMC_loop(ERGMM_MCMC_Model *model, ERGMM_MCMC_Priors *prior,
		   ERGMM_MCMC_MCMCState *cur, ERGMM_MCMC_MCMCSettings *setting, ERGMM_MCMC_ROutput *outlists);

void MBC_MCMC_store_iteration(unsigned int pos, ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par,
			      ERGMM_MCMC_MCMCSettings *setting, ERGMM_MCMC_ROutput *outlists);

#endif
