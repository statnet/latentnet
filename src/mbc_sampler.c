#include "mbc_sampler.h"
#include "ergmm_structs.h"
#include "matrix_utils.h"
#include "ergmm_probs.h"
#include "ergmm_updaters.h"


void MBC_MCMC_wrapper(int *samples_stored,
		      int *interval,

		      int *n,
		      int *d,
		      int *G,
			  
		      double *lpZList, 
		      double *lpLVList, 
		      
		      double *vZ,
		      
		      double *epsilon,
		      double *mu,
		      double *Sigma,
		      int *Ki,
		      
		      double *Sigprior, 
		      double *muSigprior, 
		      double *dirprior,
		      double *alphaprior,
		      
		      int *KiList, 
		      double *Z_pKList, 
		      double *muList, 
		      double *SigmaList){
  double **Z = vZ ? Runpack_dmatrix(vZ,*n,*d, NULL) : NULL;
  double **Z_mean_start = mu ? Runpack_dmatrix(mu,*G,*d,NULL) : NULL;
  
  /* R function enabling uniform RNG */
  GetRNGstate();
  
  
  /* Since random effects are optional (can be NULL), we have to check before
     dereferincing pointers that deal with them. */
  MBC_MCMC_init(*samples_stored, *interval,
		
		*n, *d, *G,
		
		lpZList, lpLVList,
		
		Z, 
		epsilon,Z_mean_start,Sigma,(unsigned int *) Ki,
		Sigprior? *Sigprior : 0,
		muSigprior ? *muSigprior : 0,
		dirprior ? *dirprior : 0,
		alphaprior ? *alphaprior : 0,
		KiList, Z_pKList, muList, SigmaList);
  
  PutRNGstate();
  P_free_all();
  // R's garbage collector takes care of freeing memory.
  return;
}

/* Initializes the MCMC sampler and allocates memory. */
void MBC_MCMC_init(unsigned int samples_stored, 
		   unsigned int interval, 

		   unsigned int n,
		   unsigned int d,
		   unsigned int G,
		   
		   double *lpZList,
		   double *lpLVList, 
		   
		   double **Z,

		   double *epsilon, 
		   double **Z_mean_start, 
		   double *Sigma, 
		   unsigned int *Ki,

		   double Sigprior,
		   double muSigprior,
		   double dirprior,
		   double alphaprior,

		   int *KList,
		   double *Z_pKList,
		   double *muList,
		   double *SigmaList){
  unsigned int i;


  // Packing constants into structs.
  ERGMM_MCMC_Model model = {0,
			    NULL, // Y
			    NULL, // X
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

  ERGMM_MCMC_MCMCSettings setting = {0,0,0,0,0,NULL, // deltas
				     NULL, // X_means
				     samples_stored,interval};

  ERGMM_MCMC_Priors prior = {muSigprior, // Z_mean_var
			     Sigprior, // Z_var
			     alphaprior, // a.k.a. Z_var_df (I hope)
			     NULL,
			     NULL,
			     dirprior,
			     0,0,0,0};
  
  ERGMM_MCMC_Par state = {Z, // Z
			  NULL, // coef
			  Z_mean_start, // Z_mean
			  Sigma, // Z_var
			  epsilon, // Z_pK			  
			  NULL,
			  0,
			  NULL,
			  0,
			  Ki, // Z_K
			  0, // llk
			  NULL, // lpedge
			  0, // lpZ		  
			  0, // lpLV
			  0, // lpcoef
			  0, // lpRE
			  0 // lpREV
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
				FALSE, // after_Gibbs
				NULL // update_order
  };
  
  ERGMM_MCMC_ROutput outlists = {NULL, // llk
				 lpZList,
				 NULL, // lpcoef
				 NULL, // lpRE
				 lpLVList,
				 NULL, //lpREV,
				 NULL, // Z
				 NULL, // Z_rate_move
				 NULL, // Z_rate_move_all
				 NULL, // coef
				 NULL, // coef_rate
				 muList,SigmaList,Z_pKList,
				 NULL,NULL,
				 NULL,NULL,
				 KList};

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

