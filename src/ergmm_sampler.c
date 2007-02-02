/****************************************************************************/
/*  Original Author: Susan Shortreed, susanms@stat.washington.edu           */
/*  Updated by: Jeremy Tantrum, tantrum@stat.washington.edu                 */
/*  Purpose: main functions for parameter estimation model 2                */
/*           All of this code is for an R function which is incorporated    */
/*           into the R ERGM package.                                       */
/****************************************************************************/

#include "ergmm_structs.h"
#include <R.h>
#include <R_ext/Utils.h>
#include "matrix_utils.h"
#include "ergmm_utils.h"
#include "ergmm_updaters.h"
#include "ergmm_probs.h"
#include "ergmm_families.h"
#include "wishart.h"
#include "mvnorm.h"

#include "ergmm_sampler.h"


/* Here is the wrapper function we will call from the ergm R code */
/* NOTE THE RANDOM SEED MUST BE SET IN R*/
/* Changed mostly to convert "vectorized" matrices of R to double** matrices.
 * It's tempting to also make it do the vectorizing of the result lists, but
 * that would waste (really) enormous amounts of memory.
 */

void ERGMM_MCMC_wrapper(int *samples_stored, 
			int *interval,
			   
			int *n, 
			int *p,
			int *d, 
			int *G,
			  
			int *dir,
			int *vY,
			int *family,
			int *iconsts,
			double *dconsts,

			double *vX,
			  
			double *llkList,
			double *lpZList,
			double *lpcoefList,
			double *lpREList,
			double *lpLVList,
			double *lpREVList,
			   
			double *vZ_start,

			double *epsilon,
			double *mu,
			double *Sigma,
			int *Ki,

			double *Sigprior,
			double *muSigprior, 
			double *dirprior,
			double *alphaprior,

			double *vZ_post,
			double *Z_rate_move,
			double *Z_rate_move_all,

			int *KiList,
			double *Z_pKList,
			double *muList,
			double *SigmaList,
			  
			double *coef_start,
			double *coef_mean,
			double *coef_var,
			double *Coef,
			double *B_rate, 
			  
			double *sender,
			double *receiver,
			double *sender_var,
			double *receiver_var,

			double *sender_var_prior,
			double *sender_var_prior_df,
			double *receiver_var_prior,
			double *receiver_var_prior_df,

			double *senderList,
			double *receiverList,
			double *sender_varList,
			double *receiver_varList,

			int *sociality,
			int *vobserved_ties,
			double *deltas){
  int i=0,j=0,k=0;
  double **Z_start = vZ_start ? Runpack_dmatrix(vZ_start,*n,*d, NULL) : NULL;
  double **Z_mu_start = mu ? Runpack_dmatrix(mu,*G,*d,NULL) : NULL;
  int **Y = Runpack_imatrix(vY, *n, *n, NULL);
  unsigned int **observed_ties = vobserved_ties ? (unsigned int **) Runpack_imatrix(vobserved_ties,*n,*n,NULL) : NULL;
  double ***X = (double ***) R_alloc(*p,sizeof(double**)); 

  // set up all of the covariate matrices if covariates are involed 
  // if p=0 (ie no covariates then these next two loops will do nothing)
  //

  for(i=0;i<*p;i++){
    X[i] = dmatrix((*n),(*n));
  } 
  for(k=0;k<*p;k++){
    for(i=0;i<*n;i++){
      for(j=0;j<*n;j++){
	X[k][i][j] = vX[ k*(*n)*(*n) + i*(*n) + j ];
      }
    }
  }

  /* R function enabling uniform RNG */
  GetRNGstate();
 

  /* Since random effects are optional (can be NULL), we have to check before
     dereferincing pointers that deal with them. */
  ERGMM_MCMC_init(*samples_stored, *interval,

		  *n,*p,
		  d ? *d : 0,
		  *G,
		 
		  *dir,Y,
		  family ? *family : 0,iconsts,dconsts,
		  X,

		  llkList, lpZList, lpcoefList, lpREList, lpLVList, lpREVList,
		    
		  Z_start, 
		  epsilon,Z_mu_start,Sigma,(unsigned int *)Ki,
		  Sigprior? *Sigprior : 0,
		  muSigprior ? *muSigprior : 0,
		  dirprior ? *dirprior : 0,
		  alphaprior ? *alphaprior : 0,
		  vZ_post, Z_rate_move, Z_rate_move_all, KiList, Z_pKList, muList, SigmaList,

		  coef_start,
		  Coef, B_rate,    
		  coef_mean, coef_var,

		  sender,receiver,
		  sender_var ? *sender_var : 0,
		  receiver_var ? *receiver_var : 0,
		  sender_var_prior ? *sender_var_prior : 0,
		  sender_var_prior_df ? *sender_var_prior_df : 0,
		  receiver_var_prior ? *receiver_var_prior : 0,
		  receiver_var_prior_df ? *receiver_var_prior_df : 0,
		  senderList, receiverList, 
		  sender_varList, receiver_varList,
		  *sociality,
		  observed_ties,
		  deltas[0],deltas[1],deltas[2],
		  deltas[3],deltas[4],deltas+COEF_DELTA_START);

  PutRNGstate();

  return;
}

/* Initializes the MCMC sampler and allocates memory. */
void ERGMM_MCMC_init(unsigned int samples_stored, unsigned int interval, 

		     unsigned int n, 
		     unsigned int p, unsigned int d, unsigned int G,

		     unsigned int dir, int **Y,

		     unsigned int family, int *iconsts, double *dconsts,

		     double ***X,

		     double *llkList, double *lpZList, double *lpcoefList, double *lpREList, double *lpLVList, double *lpREVList,

		     double **Z_start,
		     double *epsilon, double **Z_mu_start, double *Sigma, unsigned int *Ki,
		     double Sigprior, double muSigprior, double dirprior,
		     double alphaprior,
		     double *ZList, double *Z_rate_move, double *Z_rate_move_all, int *KList,
		     double *Z_pKList,
		     double *muList, double *SigmaList,

		     double *coef_mle,
		     double *coefList, double *coef_rate, 
		     double *coef_mean, double *coef_var, 

		     double *sender, double *receiver,
		     double sender_var, double receiver_var,
		     double sender_var_prior, double sender_var_prior_df,
		     double receiver_var_prior, double receiver_var_prior_df,
		     double *senderList, double *receiverList,
		     double *sender_varList, double *receiver_varList,
		     unsigned int sociality,
		     unsigned int **observed_ties,

		     double Z_delta, double Z_tr_delta, double Z_scl_delta,
		     double RE_delta, double RE_shift_delta,
		     double *coef_delta)
{
  unsigned int i,j,k,n_observed;
  double X_sum;

  // Packing constants into structs.
  ERGMM_MCMC_Model model = {dir,
			    Y, // Y
			    X, // X
			    (unsigned int **)observed_ties,
			    NULL,
			    0,
			    iconsts,
			    dconsts,
			    n, // verts
			    d, // latent
			    p, // coef
			    G, // clusters
			    sociality};

  // Set the p(Y_ij|eta_ij) function and compute the term constant in eta.
  switch(family){
  case 0:
    model.lp_edge=ERGMM_MCMC_lp_edge_Bernoulli_logit;
    ERGMM_MCMC_set_lp_Yconst_Bernoulli_logit(&model);
    break;
  case 1:
    model.lp_edge=ERGMM_MCMC_lp_edge_binomial_logit;
    ERGMM_MCMC_set_lp_Yconst_binomial_logit(&model);
    break;
  case 2:
    model.lp_edge=ERGMM_MCMC_lp_edge_Poisson_log;
    ERGMM_MCMC_set_lp_Yconst_Poisson_log(&model);
    break;
  }

  ERGMM_MCMC_MCMCSettings setting = {Z_delta,Z_tr_delta,Z_scl_delta,
				     RE_delta,RE_shift_delta,coef_delta,
				     dvector(model.coef), // X_means
				     samples_stored,interval};

  for(k=0;k<model.coef;k++){
    X_sum=0;
    n_observed=0;
    for(i=0;i<model.verts;i++){
      if(model.dir)
	for(j=0;j<model.verts;j++){
	  n_observed+=IS_OBSERVABLE(model.observed_ties,i,j) ? 1:0;
	  X_sum+=model.X[k][i][j]*IS_OBSERVABLE(model.observed_ties,i,j);
	}
      else
	for(j=0;j<i;j++){
	  n_observed+=IS_OBSERVABLE(model.observed_ties,i,j) ? 1:0;
	  X_sum+=model.X[k][i][j]*IS_OBSERVABLE(model.observed_ties,i,j);
	}
    }
    setting.X_means[k]=X_sum/n_observed;
  }

  ERGMM_MCMC_Priors prior = {muSigprior, // Z_mu_var
			     Sigprior, // Z_var
			     alphaprior, // a.k.a. Z_var_df (I hope)
			     coef_mean,
			     coef_var,
			     dirprior,
			     sender_var_prior,
			     sender_var_prior_df,
			     receiver_var_prior,
			     receiver_var_prior_df};
  
  ERGMM_MCMC_Par state = {Z_start, // Z
			  coef_mle, // coef
			  Z_mu_start, // Z_mu
			  Sigma, // Z_var
			  epsilon, // Z_pK			  
			  sender,
			  sender_var,
			  model.sociality?sender:receiver,
			  receiver_var,
			  Ki, // Z_K
			  0, // llk
			  dmatrix(model.verts,model.verts), // lpedge
			  0, // lpZ		  
			  0, // lpLV
			  0, // lpcoef
			  0, // lpRE
			  0 // lpREV
  };

  ERGMM_MCMC_Par prop = {model.latent ? dmatrix(model.verts,model.latent):NULL, // Z
			 model.coef ? dvector(model.coef):NULL, // coef
			 model.clusters ? dmatrix(model.clusters,model.latent):NULL, // Z_mu
			 model.latent ? dvector(model.clusters?model.clusters:1):NULL, // Z_var
			 model.clusters ? dvector(model.clusters):NULL, // Z_pK
			 sender ? dvector(model.verts):NULL, // sender
			 0, // sender_var
			 receiver && !model.sociality ? dvector(model.verts):NULL, // receiver
			 0,
			 state.Z_K, // prop.Z_K === state.Z_K
			 0, // llk
			 dmatrix(model.verts,model.verts), // lpedge
			 0, // lpZ
			 0, // lpLV
			 0, // lpcoef
			 0, // lpRE
			 0 // lpREV
  };
  if(model.sociality) prop.receiver=prop.sender;

  ERGMM_MCMC_MCMCState start = {&state,
				&prop,
				model.clusters ? dmatrix(model.clusters,model.latent) : NULL, // Z_bar
				model.latent ? dvector(model.latent) : NULL, // tr_by
				model.clusters ? dvector(model.clusters): NULL, // pK
				model.clusters ? (unsigned int *) ivector(model.clusters) : NULL, // n
				PROP_NONE, // prop_Z
				PROP_NONE, // prop_RE
				PROP_NONE, // prop_coef
				PROP_NONE, // prop_LV
				PROP_NONE, // prop_REV
				FALSE, // after_Gibbs
				(model.latent || sender || receiver) ? (unsigned int *) ivector(model.verts) : NULL // update_order
  };
  
  ERGMM_MCMC_ROutput outlists = {llkList, lpZList, lpcoefList, lpREList, lpLVList, lpREVList,
				 ZList, Z_rate_move, Z_rate_move_all,
				 coefList,coef_rate,
				 muList,SigmaList,Z_pKList,
				 senderList,sender_varList,
				 receiverList,receiver_varList,
				 KList};

  if(model.clusters>0)
    for(i=0;i<model.verts;i++)
      start.n[state.Z_K[i] - 1]++;

  // Initialize the log-probabilities.
  state.llk = ERGMM_MCMC_lp_Y(&model, &state,0);
  if(model.latent) ERGMM_MCMC_logp_Z(&model, &state);
  if(state.sender || state.receiver) ERGMM_MCMC_logp_RE(&model, &state);
  if(state.coef) ERGMM_MCMC_logp_coef(&model, &state, &prior);
  if(model.latent) ERGMM_MCMC_logp_LV(&model, &state, &prior);
  if(state.sender || state.receiver) ERGMM_MCMC_logp_REV(&model, &state, &prior);
  copy_MCMC_Par(&model,&state,&prop);
  ERGMM_MCMC_store_iteration(0,&model,&state,&setting,&outlists);
  ERGMM_MCMC_store_iteration(1,&model,&state,&setting,&outlists);

  ERGMM_MCMC_loop(&model, &prior, &start, &setting, &outlists);


  // R's garbage collector takes care of freeing memory.
}


void ERGMM_MCMC_loop(ERGMM_MCMC_Model *model, ERGMM_MCMC_Priors *prior,
		     ERGMM_MCMC_MCMCState *cur, ERGMM_MCMC_MCMCSettings *setting, ERGMM_MCMC_ROutput *outlists)
{  
  unsigned int n_accept_z=0, n_accept_transl_z=0, n_accept_b=0, pos=0;
  unsigned int iter, total_iters = setting->sample_size*setting->interval;

  /* Note that indexing here starts with 1.
     It can be thought of as follows:
     At the end of the updates, we have made iter complete MCMC updates. */
  for(iter=1;iter<=total_iters;iter++){

    R_CheckUserInterrupt(); // So that CTRL-C can interrupt the run.

    n_accept_z += ERGMM_MCMC_Z_RE_up(model, prior, cur, setting); 

    if(model->latent){

      //n_accept_transl_z += translate_Z(model,prior,cur,setting);

      // Update cluster parameters (they are separated from data by Z, so plain Gibbs).
      // Note that they are also updated in coef_up_scl_tr_Z_shift_RE.
      if(model->clusters>0)
	ERGMM_MCMC_CV_up(model,prior,cur);
      else
	ERGMM_MCMC_LV_up(model,prior,cur);
    }

    /* Update coef given this new value of Z and conditioned on everything else.
       Also propose to scale Z and shift random effects.
    */
    
    if( ERGMM_MCMC_coef_up_scl_tr_Z_shift_RE(model,prior,cur,setting) ){
      n_accept_b++;
    }  

    if(cur->state->sender || cur->state->receiver){
      ERGMM_MCMC_REV_up(model,prior,cur);
    }

    // If we have a new MLE (actually, the highest likelihood encountered to date), store it.
    if( cur->state->llk > GET_DEFAULT(outlists->llk,0,0) ) ERGMM_MCMC_store_iteration(0,model,cur->state,setting,outlists);

    // If we have a new posterior mode (actually, the highest joint density of all variables but K observed to date), store it.
    if( cur->state->llk + cur->state->lpZ + cur->state->lpLV + 
	cur->state->lpcoef + cur->state->lpRE + cur->state->lpREV >
	GET_DEFAULT(outlists->llk,1,0) + GET_DEFAULT(outlists->lpZ,1,0) + GET_DEFAULT(outlists->lpLV,1,0) + 
	GET_DEFAULT(outlists->lpcoef,1,0) + GET_DEFAULT(outlists->lpRE,1,0) + GET_DEFAULT(outlists->lpREV,1,0) )
      ERGMM_MCMC_store_iteration(1,model,cur->state,setting,outlists);

    /* every interval save the results */
    if((iter % setting->interval) == 0){
      pos = (iter/setting->interval-1)+ERGMM_OUTLISTS_RESERVE;

      // Store current iteration.
      ERGMM_MCMC_store_iteration(pos, model, cur->state, setting, outlists);

      // Acceptance rates.
      outlists->coef_rate[pos] = (double) ((double)n_accept_b)/((double)setting->interval);
      if(model->latent){
	outlists->Z_rate_move[pos] = (double) ((double)n_accept_z)/((double)setting->interval*model->verts); 
	outlists->Z_rate_move_all[pos] = (double) ((double)n_accept_transl_z)/((double)setting->interval); 
      }

      n_accept_z=0; 
      n_accept_b=0; 
      n_accept_transl_z=0;
    }

  } // end main MCMC loop

  return;
  
}

void ERGMM_MCMC_store_iteration(unsigned int pos, ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par,
		     ERGMM_MCMC_MCMCSettings *setting, ERGMM_MCMC_ROutput *outlists){
  // Log-likelihood and other probabilities:
  outlists->llk[pos] = par->llk; 
  if(outlists->lpZ)
    outlists->lpZ[pos] = par->lpZ;
  if(outlists->lpcoef)
    outlists->lpcoef[pos] = par->lpcoef;
  if(outlists->lpRE)
    outlists->lpRE[pos] = par->lpRE;
  if(outlists->lpLV)
    outlists->lpLV[pos] = par->lpLV;
  if(outlists->lpREV)
    outlists->lpREV[pos] = par->lpREV;

  // Covariate coefficients.
  Rpack_dvector(par->coef,model->coef,outlists->coef+pos,setting->sample_size+ERGMM_OUTLISTS_RESERVE);


  if(model->latent){
    // Latent positions.
    Rpack_dmatrix(par->Z,model->verts,model->latent,outlists->Z+pos,setting->sample_size+ERGMM_OUTLISTS_RESERVE);

      
    // Cluster-related.
    if(model->clusters>0){      
      // Cluster assignments.
      Rpack_ivector((int *)par->Z_K,model->verts,outlists->Z_K+pos,setting->sample_size+ERGMM_OUTLISTS_RESERVE);
	  
      // Cluster means.
      Rpack_dmatrix(par->Z_mu,model->clusters,model->latent,outlists->Z_mu+pos,setting->sample_size+ERGMM_OUTLISTS_RESERVE);
	  
      // Intracluster variances.
      Rpack_dvector(par->Z_var,model->clusters,outlists->Z_var+pos,setting->sample_size+ERGMM_OUTLISTS_RESERVE);

      // Cluster probabilities.
      Rpack_dvector(par->Z_pK,model->clusters,outlists->Z_pK+pos,setting->sample_size+ERGMM_OUTLISTS_RESERVE);

    }
    else
      outlists->Z_var[pos]=par->Z_var[0];
  }

  // Sender effects.
  if(par->sender){
    Rpack_dvector(par->sender,model->verts,outlists->sender+pos,setting->sample_size+ERGMM_OUTLISTS_RESERVE);
    outlists->sender_var[pos] = par->sender_var;
  }

  // Receiver effects.
  if(par->receiver && !model->sociality){
    Rpack_dvector(par->receiver,model->verts,outlists->receiver+pos,setting->sample_size+ERGMM_OUTLISTS_RESERVE);
    outlists->receiver_var[pos] = par->receiver_var;
  }      
      
}

