/*  File src/mbc_sampler.c in package latentnet, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2018 Statnet Commons
 */
/*********************************************/
/* Bayesian model-based clustering routines. */
/*********************************************/

#include "mbc_sampler.h"
#include "ergmm_structs.h"
#include "matrix_utils.h"
#include "ergmm_probs.h"
#include "ergmm_updaters.h"


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
		      double *Z_var_mcmc){
  double **Z = Runpack_dmatrix(vZ,*n,*d, NULL);
  double **Z_mean_start = Runpack_dmatrix(Z_mean,*G,*d,NULL);
  
  /* R function enabling uniform RNG */
  GetRNGstate();
  
  
  MBC_MCMC_init(*sample_size, *interval,
		
		*n, *d, *G,
		
		lpZ_mcmc, lpLV_mcmc,
		
		Z, 
		Z_pK,Z_mean_start,Z_var,(unsigned int *) Z_K,
		*Z_var_prior,
		*Z_mean_prior_var,
		*Z_K_prior,
		*Z_var_df,
		Z_K_mcmc, Z_pK_mcmc, Z_mean_mcmc, Z_var_mcmc);
  
  PutRNGstate();
  P_free_all();
  // R's garbage collector takes care of freeing memory.
  return;
}

/* Initializes the MCMC sampler and allocates memory. */
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
		   double *Z_var_mcmc){
  unsigned int i;


  // Packing constants into structs.
  ERGMM_MCMC_Model model = {0,
			    NULL, // iY
			    NULL, // dY
			    NULL, // X
			    NULL,
			    NULL,
			    NULL,
			    0,
			    NULL,
			    NULL,
			    n, // verts
			    d, // latent
			    0, // coef
			    G // clusters
  };

  ERGMM_MCMC_MCMCSettings setting = {0,0,NULL,NULL,NULL,0,0,0, // deltas
				     sample_size,interval,
				     FALSE // accept_all
  };

  ERGMM_MCMC_Priors prior = {Z_mean_prior_var, // Z_mean_var
			     Z_var_prior, // Z_var
			     Z_var_df, // a.k.a. Z_var_df (I hope)
			     NULL,
			     NULL,
			     Z_K_prior,
			     0,0,0,0};
  
  ERGMM_MCMC_Par state = {Z, // Z
			  NULL, // coef
			  Z_mean_start, // Z_mean
			  Z_var, // Z_var
			  Z_pK, // Z_pK			  
			  NULL,
			  0,
			  NULL,
			  0,
			  0, // dispersion
			  Z_K, // Z_K
			  0, // llk
			  NULL, // lpedge
			  0, // lpZ		  
			  0, // lpLV
			  0, // lpcoef
			  0, // lpRE
			  0, // lpREV
			  0  // lpdispersion
  };

  ERGMM_MCMC_MCMCState start = {&state,
				NULL,
				model.clusters ? dmatrix(model.clusters,model.latent) : NULL, // Z_bar
				NULL, // tr_by
				model.clusters ? dvector(model.clusters): NULL, // pK
				model.clusters ? (unsigned int *) ivector(model.clusters) : NULL, // n
				PROP_NONE, // prop_Z
				PROP_NONE, // prop_RE
				PROP_NONE, // prop_coef
				PROP_NONE, // prop_LV
				PROP_NONE, // prop_REV
				PROP_NONE, // prop_dispersion
				FALSE, // after_Gibbs
				NULL // update_order
  };
  
  ERGMM_MCMC_ROutput outlists = {NULL, // llk
				 lpZ_mcmc,
				 NULL, // lpcoef
				 NULL, // lpRE
				 lpLV_mcmc,
				 NULL, //lpREV,
				 NULL, //lpdispersion
				 NULL, // Z
				 NULL, // Z_rate_move
				 NULL, // coef
				 NULL, // coef_rate
				 Z_mean_mcmc,Z_var_mcmc,Z_pK_mcmc,
				 NULL,NULL,
				 NULL,NULL,
				 NULL, // dispersion_mcmc
				 K_mcmc};

  if(model.clusters>0)
    for(i=0;i<model.verts;i++)
      start.n[state.Z_K[i] - 1]++;

  // Initialize the log-probabilities.
  ERGMM_MCMC_logp_Z(&model, &state);
  MBC_MCMC_store_iteration(0,&model,&state,&setting,&outlists);
  MBC_MCMC_store_iteration(1,&model,&state,&setting,&outlists);
  
  MBC_MCMC_loop(&model, &prior, &start, &setting, &outlists);

}


void MBC_MCMC_loop(ERGMM_MCMC_Model *model, ERGMM_MCMC_Priors *prior,
		   ERGMM_MCMC_MCMCState *cur, ERGMM_MCMC_MCMCSettings *setting, ERGMM_MCMC_ROutput *outlists)
{  
  unsigned int pos=0;
  unsigned int iter, total_iters = setting->sample_size*setting->interval;

  //  Rprintf("Started MCMC loop.\n");
  /* Note that indexing here starts with 1.
     It can be thought of as follows:
     At the end of the updates, we have made iter complete MCMC updates. */
  for(iter=1;iter<=total_iters;iter++){

    R_CheckUserInterrupt(); // So that CTRL-C can interrupt the run.

    ERGMM_MCMC_CV_up(model,prior,cur);
    ERGMM_MCMC_logp_Z(model, cur->state);
    
    // If we have a new MLE (actually, the highest likelihood encountered to date), store it.
    if( cur->state->lpZ > outlists->lpZ[0] ) MBC_MCMC_store_iteration(0,model,cur->state,setting,outlists);
    if( cur->state->lpZ + cur->state->lpLV > outlists->lpZ[1] + outlists->lpLV[1] ) MBC_MCMC_store_iteration(1,model,cur->state,setting,outlists);

    /* every interval save the results */
    if((iter % setting->interval) == 0){
      pos = (iter/setting->interval-1)+MBC_OUTLISTS_RESERVE;

      // Store current iteration.
      MBC_MCMC_store_iteration(pos, model, cur->state, setting, outlists);

    }

  } // end main MCMC loop

  return;
  
}

void MBC_MCMC_store_iteration(unsigned int pos, ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par,
			      ERGMM_MCMC_MCMCSettings *setting, ERGMM_MCMC_ROutput *outlists){
  // Log-likelihood a.k.a. lpZ.
  outlists->lpZ[pos] = par->lpZ; 
  outlists->lpLV[pos] = par->lpLV; 
  
  // Cluster-related.
  // Cluster assignments.
  Rpack_ivectors((int *)par->Z_K,model->verts,outlists->Z_K+pos,setting->sample_size+MBC_OUTLISTS_RESERVE);
  
  // Cluster means.
  Rpack_dmatrixs(par->Z_mean,model->clusters,model->latent,outlists->Z_mean+pos,setting->sample_size+MBC_OUTLISTS_RESERVE);
  
  // Intracluster variances.
  Rpack_dvectors(par->Z_var,model->clusters,outlists->Z_var+pos,setting->sample_size+MBC_OUTLISTS_RESERVE);
  
  // Cluster probabilities.
  Rpack_dvectors(par->Z_pK,model->clusters,outlists->Z_pK+pos,setting->sample_size+MBC_OUTLISTS_RESERVE);
  
}

