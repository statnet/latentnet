/*  File src/ergmm_sampler.c in package latentnet, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2020 Statnet Commons
 */
/****************************************************************************/
/* Top-level functions for sampling from an ERGMM's posterior distribution. */
/****************************************************************************/

#include "ergmm_structs.h"
#include <R.h>
#include <R_ext/Utils.h>
#include "matrix_utils.h"
#include "ergmm_utils.h"
#include "ergmm_updaters.h"
#include "ergmm_probs.h"
#include "ergmm_families.h"
#include "ergmm_latent_effects.h"
#include "wishart.h"
#include "mvnorm.h"

#include "ergmm_sampler.h"


/* Here is the wrapper function we will call from the ergm R code */
/* NOTE THE RANDOM SEED MUST BE SET IN R*/
/* Changed mostly to convert "vectorized" matrices of R to double** matrices.
 * It's tempting to also make it do the vectorizing of the result lists, but
 * that would waste (really) enormous amounts of memory.
 */

void ERGMM_MCMC_wrapper(int *sample_size, 
			int *interval,
			   
			int *n, 
			int *p,
			int *d, 
			int *G,
			int *latent_eff,
			int *family,
			int *res,
			  
			int *dir,
			int *viY,
			double *vdY,
			int *iconsts,
			double *dconsts,

			double *vX,
			  
			double *llk_mcmc,
			double *lpZ_mcmc,
			double *lpcoef_mcmc,
			double *lpRE_mcmc,
			double *lpLV_mcmc,
			double *lpREV_mcmc,
			double *lpdispersion_mcmc,
			   
			double *vZ_start,

			double *Z_pK_start,
			double *vZ_mean_start,
			double *Z_var_start,
			int *Z_K_start,

			double *Z_var_prior,
			double *Z_mean_prior_var, 
			double *Z_pK_prior,
			double *Z_var_prior_df,

			double *Z_mcmc,
			double *Z_rate_move,

			int *Z_K_mcmc,
			double *Z_pK_mcmc,
			double *Z_mean_mcmc,
			double *Z_var_mcmc,
			  
			double *coef_start,
			double *coef_prior_mean,
			double *coef_var,
			double *coef_mcmc,
			double *coef_rate, 
			  
			double *sender_start,
			double *receiver_start,
			double *sender_var_start,
			double *receiver_var_start,

			double *sender_var_prior,
			double *sender_var_prior_df,
			double *receiver_var_prior,
			double *receiver_var_prior_df,

			double *sender_mcmc,
			double *receiver_mcmc,
			double *sender_var_mcmc,
			double *receiver_var_mcmc,

			double *dispersion_start,
			double *dispersion_prior,
			double *dispersion_prior_df,
			double *dispersion_mcmc,

			int *vobserved_ties,

			double *deltas,
			double *vcoef_eff_sender,
			int *coef_eff_sender_size,
			double *vcoef_eff_receiver,
			int *coef_eff_receiver_size,

			int *accept_all){

  // This was added because newer versions of R no longer pass a 0-length vector as NULL, so we have to do it here.
  if(*p==0){
    vX=NULL; lpcoef_mcmc=NULL; coef_start=NULL; coef_prior_mean=NULL; coef_var=NULL; coef_mcmc=NULL; coef_rate=NULL; vcoef_eff_sender=NULL; coef_eff_sender_size=NULL; vcoef_eff_receiver=NULL; coef_eff_receiver_size=NULL;
  }
  if(*G==0){
    Z_pK_start=NULL; vZ_mean_start=NULL; Z_K_start=NULL; Z_mean_prior_var=NULL; Z_pK_prior=NULL; Z_K_mcmc=NULL; Z_pK_mcmc=NULL; Z_mean_mcmc=NULL;
  }
  if(*d==0){
    lpZ_mcmc=NULL; lpLV_mcmc=NULL; vZ_start=NULL; Z_pK_start=NULL; vZ_mean_start=NULL; Z_var_start=NULL; Z_K_start=NULL; Z_var_prior=NULL; Z_mean_prior_var=NULL; Z_pK_prior=NULL; Z_var_prior_df=NULL; Z_mcmc=NULL; Z_K_mcmc=NULL; Z_pK_mcmc=NULL; Z_mean_mcmc=NULL; Z_var_mcmc=NULL;
  }
  if(res[0]==0&&res[2]==0){
    sender_start=NULL; sender_var_start=NULL; sender_var_prior=NULL; sender_var_prior_df=NULL; sender_mcmc=NULL; sender_var_mcmc=NULL; vcoef_eff_sender=NULL; coef_eff_sender_size=NULL;
  }
  if(res[1]==0){
    receiver_start=NULL; receiver_var_start=NULL; receiver_var_prior=NULL; receiver_var_prior_df=NULL; receiver_mcmc=NULL; receiver_var_mcmc=NULL; vcoef_eff_receiver=NULL; coef_eff_receiver_size=NULL;
  }
  if(res[0]==0&&res[1]==0&&res[2]==0){
    lpRE_mcmc=NULL; lpREV_mcmc=NULL;
    if(*d==0) Z_rate_move=NULL;
  }
  if(res[3]==0){
    lpdispersion_mcmc=NULL;
    dispersion_mcmc=NULL;
  }


  int i=0,j=0,k=0;
  double **Z_start = vZ_start ? Runpack_dmatrix(vZ_start,*n,*d, NULL) : NULL;
  double **Z_mean_start = vZ_mean_start ? Runpack_dmatrix(vZ_mean_start,*G,*d,NULL) : NULL;
  int **iY = Runpack_imatrix(viY, *n, *n, NULL);
  double **dY = Runpack_dmatrix(vdY, *n, *n, NULL);
  unsigned int **observed_ties = *vobserved_ties>=0 ? (unsigned int **) Runpack_imatrix(vobserved_ties,*n,*n,NULL) : NULL;
  double ***X = d3array(*p,*n,*n);

  /* The joint proposal coefficient matrix is square with side
     + covariate coefficients  : p
     + latent space            : 1
     + sender                  : ~p
     + receiver (no sociality) : ~p
     + dispersion              : 1
  */

  unsigned int group_prop_size = *p + (*d ? 1 : 0) + (vcoef_eff_sender?*coef_eff_sender_size:0)+(vcoef_eff_receiver?*coef_eff_receiver_size:0) + (lpdispersion_mcmc?1:0);
  double **group_deltas = Runpack_dmatrix(deltas+GROUP_DELTAS_START, group_prop_size, group_prop_size, NULL);
  double **coef_eff_sender = vcoef_eff_sender? Runpack_dmatrix(vcoef_eff_sender, *coef_eff_sender_size, *n, NULL) : NULL;
  double **coef_eff_receiver = vcoef_eff_receiver? Runpack_dmatrix(vcoef_eff_receiver, *coef_eff_receiver_size, *n, NULL) : NULL;


  // set up all of the covariate matrices if covariates are involed 
  // if p=0 (ie no covariates then these next two loops will do nothing)
  //

  for(k=0;k<*p;k++){
    for(i=0;i<*n;i++){
      for(j=0;j<*n;j++){
	X[k][i][j] = vX[ k*(*n)*(*n) + j*(*n) + i ];
      }
    }
  }

  /* R function enabling uniform RNG */
  GetRNGstate();
 
  ERGMM_MCMC_init(*sample_size, *interval,

		  *n,*p,
		  d ? *d : 0,
		  *G,
		 
		  *dir,iY,dY,
		  *family-1,
		  iconsts,dconsts,
		  *latent_eff ? *latent_eff-1 : 1000,

		  X,

		  llk_mcmc, lpZ_mcmc, lpcoef_mcmc, lpRE_mcmc, lpLV_mcmc, lpREV_mcmc, lpdispersion_mcmc,
		    
		  Z_start, 
		  Z_pK_start,Z_mean_start,Z_var_start,(unsigned int *)Z_K_start,
		  Z_var_prior? *Z_var_prior : 0,
		  Z_mean_prior_var ? *Z_mean_prior_var : 0,
		  Z_pK_prior ? *Z_pK_prior : 0,
		  Z_var_prior_df ? *Z_var_prior_df : 0,
		  Z_mcmc, Z_rate_move, Z_K_mcmc, Z_pK_mcmc, Z_mean_mcmc, Z_var_mcmc,

		  coef_start,
		  coef_mcmc, coef_rate,    
		  coef_prior_mean, coef_var,

		  sender_start,receiver_start,
		  sender_var_start ? *sender_var_start : 0,
		  receiver_var_start ? *receiver_var_start : 0,
		  sender_var_prior ? *sender_var_prior : 0,
		  sender_var_prior_df ? *sender_var_prior_df : 0,
		  receiver_var_prior ? *receiver_var_prior : 0,
		  receiver_var_prior_df ? *receiver_var_prior_df : 0,
		  sender_mcmc, receiver_mcmc, 
		  sender_var_mcmc, receiver_var_mcmc,

		  dispersion_mcmc? *dispersion_start:0,
		  dispersion_mcmc? *dispersion_prior:0,
		  dispersion_mcmc? *dispersion_prior_df:0,
		  dispersion_mcmc,

		  res[2],
		  observed_ties,
		  deltas[0],deltas[1],group_deltas,group_prop_size,
		  coef_eff_sender,coef_eff_sender_size?*coef_eff_sender_size:0,
		  coef_eff_receiver,coef_eff_receiver_size?*coef_eff_receiver_size:0,
		  *accept_all);

  PutRNGstate();
  P_free_all();
  return;
}

/* Initializes the MCMC sampler and allocates memory. */
void ERGMM_MCMC_init(unsigned int sample_size, unsigned int interval, 

		     unsigned int n,
		     unsigned int p, unsigned int d, unsigned int G,

		     unsigned int dir, int **iY, double **dY,

		     unsigned int family, 
		     int *iconsts, double *dconsts,
		     unsigned int latent_eff,

		     double ***X,

		     double *llk_mcmc, double *lpZ_mcmc, double *lpcoef_mcmc, double *lpRE_mcmc, double *lpLV_mcmc, double *lpREV_mcmc, double *lpdispersion_mcmc,

		     double **Z_start,
		     double *Z_pK_start, double **Z_mean_start, double *Z_var_start, unsigned int *Z_K_start,
		     double Z_var_prior, double Z_mean_prior_var, double Z_pK_prior,
		     double Z_var_prior_df,
		     double *Z_mcmc, double *Z_rate_move, int *K_mcmc,
		     double *Z_pK_mcmc,
		     double *Z_mean_mcmc, double *Z_var_mcmc,

		     double *coef_start,
		     double *coef_mcmc, double *coef_rate, 
		     double *coef_prior_mean, double *coef_var, 

		     double *sender_start, double *receiver_start,
		     double sender_var_start, double receiver_var_start,
		     double sender_var_prior, double sender_var_prior_df,
		     double receiver_var_prior, double receiver_var_prior_df,
		     double *sender_mcmc, double *receiver_mcmc,
		     double *sender_var_mcmc, double *receiver_var_mcmc,

		     double dispersion_start, double dispersion_prior,
		     double dispersion_prior_df, double *dispersion_mcmc,

		     unsigned int sociality,
		     unsigned int **observed_ties,
		     
		     double Z_delta,
		     double RE_delta,
		     double **group_deltas,
		     unsigned int group_prop_size,
		     double **coef_eff_sender,
		     unsigned int coef_eff_sender_size,
		     double **coef_eff_receiver,
		     unsigned int coef_eff_receiver_size,

		     unsigned int accept_all)
{
  unsigned int i;

  // Packing constants into structs.
  ERGMM_MCMC_Model model = {dir,
			    iY, // iY
			    dY, // dY
			    X, // X
			    (unsigned int **)observed_ties,
			    ERGMM_MCMC_lp_edge[family],
			    ERGMM_MCMC_E_edge[family],
			    0,
			    iconsts,
			    dconsts,
			    n, // verts
			    d, // latent
			    p, // coef
			    G, // clusters
			    sociality,
			    dispersion_mcmc!=NULL,
			    d ? ERGMM_MCMC_latent_eff[latent_eff] : NULL
  };
  ERGMM_MCMC_set_lp_Yconst[family](&model);

  ERGMM_MCMC_MCMCSettings setting = {Z_delta,
				     RE_delta,
				     group_deltas,
				     coef_eff_sender,
				     coef_eff_receiver,
				     group_prop_size,
				     coef_eff_sender_size,
				     coef_eff_receiver_size,
				     sample_size,interval,
				     accept_all
  };

  ERGMM_MCMC_Priors prior = {Z_mean_prior_var, // Z_mean_var
			     Z_var_prior, // Z_var
			     Z_var_prior_df, // a.k.a. Z_var_df (I hope)
			     coef_prior_mean,
			     coef_var,
			     Z_pK_prior,
			     sender_var_prior,
			     sender_var_prior_df,
			     receiver_var_prior,
			     receiver_var_prior_df,
			     dispersion_prior,
			     dispersion_prior_df};
  
  ERGMM_MCMC_Par state = {Z_start, // Z
			  coef_start, // coef
			  Z_mean_start, // Z_mean
			  Z_var_start, // Z_var
			  Z_pK_start, // Z_pK			  
			  sender_start,
			  sender_var_start,
			  model.sociality?sender_start:receiver_start,
			  receiver_var_start,
			  dispersion_start,
			  Z_K_start, // Z_K
			  0, // llk
			  dmatrix(model.verts,model.verts), // lpedge
			  0, // lpZ		  
			  0, // lpLV
			  0, // lpcoef
			  0, // lpRE
			  0, // lpREV
			  0, // lpdispersion
  };

  ERGMM_MCMC_Par prop = {model.latent ? dmatrix(model.verts,model.latent):NULL, // Z
			 model.coef ? dvector(model.coef):NULL, // coef
			 model.clusters ? dmatrix(model.clusters,model.latent):NULL, // Z_mean
			 model.latent ? dvector(model.clusters?model.clusters:1):NULL, // Z_var
			 model.clusters ? dvector(model.clusters):NULL, // Z_pK
			 sender_start ? dvector(model.verts):NULL, // sender
			 0, // sender_var
			 receiver_start && !model.sociality ? dvector(model.verts):NULL, // receiver
			 0, // receiver_var
			 0, // dispersion
			 state.Z_K, // prop.Z_K === state.Z_K
			 0, // llk
			 dmatrix(model.verts,model.verts), // lpedge
			 0, // lpZ
			 0, // lpLV
			 0, // lpcoef
			 0, // lpRE
			 0, // lpREV
			 0 // lpdispersion
  };
  if(model.sociality) prop.receiver=prop.sender;

  ERGMM_MCMC_MCMCState start = {&state,
				&prop,
				model.clusters ? dmatrix(model.clusters,model.latent) : NULL, // Z_bar
				setting.group_prop_size ? dvector(setting.group_prop_size) : NULL, // deltas
				model.clusters ? dvector(model.clusters): NULL, // pK
				model.clusters ? (unsigned int *) ivector(model.clusters) : NULL, // n
				PROP_NONE, // prop_Z
				PROP_NONE, // prop_RE
				PROP_NONE, // prop_coef
				PROP_NONE, // prop_LV
				PROP_NONE, // prop_REV
				PROP_NONE, // prop_dispersion
				FALSE, // after_Gibbs
				(model.latent || sender_start || receiver_start) ? (unsigned int *) ivector(model.verts) : NULL // update_order
  };
  
  ERGMM_MCMC_ROutput outlists = {llk_mcmc, lpZ_mcmc, lpcoef_mcmc, lpRE_mcmc, lpLV_mcmc, lpREV_mcmc, lpdispersion_mcmc,
				 Z_mcmc, Z_rate_move,
				 coef_mcmc,coef_rate,
				 Z_mean_mcmc,Z_var_mcmc,Z_pK_mcmc,
				 sender_mcmc,sender_var_mcmc,
				 receiver_mcmc,receiver_var_mcmc,
				 dispersion_mcmc,
				 K_mcmc};

  if(model.clusters>0)
    for(i=0;i<model.verts;i++)
      start.n[state.Z_K[i] - 1]++;

  // Initialize the log-probabilities.
  state.llk = ERGMM_MCMC_lp_Y(&model, &state, TRUE);
  if(model.latent) ERGMM_MCMC_logp_Z(&model, &state);
  if(state.sender || state.receiver) ERGMM_MCMC_logp_RE(&model, &state);
  if(state.coef) ERGMM_MCMC_logp_coef(&model, &state, &prior);
  if(model.latent) ERGMM_MCMC_logp_LV(&model, &state, &prior);
  if(state.sender || state.receiver) ERGMM_MCMC_logp_REV(&model, &state, &prior);
  if(model.dispersion) ERGMM_MCMC_logp_dispersion(&model, &state, &prior);
  copy_MCMC_Par(&model,&state,&prop);
  ERGMM_MCMC_store_iteration(0,&model,&state,&setting,&outlists);
  ERGMM_MCMC_store_iteration(1,&model,&state,&setting,&outlists);

  ERGMM_MCMC_loop(&model, &prior, &start, &setting, &outlists);


  // R's garbage collector takes care of freeing memory.
}


void ERGMM_MCMC_loop(ERGMM_MCMC_Model *model, ERGMM_MCMC_Priors *prior,
		     ERGMM_MCMC_MCMCState *cur, ERGMM_MCMC_MCMCSettings *setting, ERGMM_MCMC_ROutput *outlists)
{  
  unsigned int n_accept_z=0, n_accept_b=0, pos=0;
  unsigned int iter, total_iters = setting->sample_size*setting->interval;

  /* Note that indexing here starts with 1.
     It can be thought of as follows:
     At the end of the updates, we have made iter complete MCMC updates. */
  for(iter=1;iter<=total_iters;iter++){

    R_CheckUserInterrupt(); // So that CTRL-C can interrupt the run.
    if(model->latent || cur->state->sender || cur->state->receiver)
      n_accept_z += ERGMM_MCMC_Z_RE_up(model, prior, cur, setting);

    if(model->latent){
      // Update cluster parameters (they are separated from data by Z, so full conditional sampling).
      // Note that they are also updated in coef_up_scl_tr_Z_shift_RE.
      if(model->clusters>0)
	ERGMM_MCMC_CV_up(model,prior,cur);
      else
	ERGMM_MCMC_LV_up(model,prior,cur);
    }

    /* Update coef given this new value of Z and conditioned on everything else.
       Also propose to scale Z and shift random effects.
    */
    if( ERGMM_MCMC_coef_up_scl_Z_shift_RE(model,prior,cur,setting) ){
      n_accept_b++;
    }

    if(cur->state->sender || cur->state->receiver){
      ERGMM_MCMC_REV_up(model,prior,cur);
    }

    // If we have a new MLE (actually, the highest likelihood encountered to date), store it.
    if( cur->state->llk > GET_DEFAULT(outlists->llk,0,0) ) ERGMM_MCMC_store_iteration(0,model,cur->state,setting,outlists);

    // If we have a new posterior mode (actually, the highest joint density of all variables but K observed to date), store it.
    if( cur->state->llk + cur->state->lpZ + cur->state->lpLV + 
	cur->state->lpcoef + cur->state->lpRE + cur->state->lpREV + cur->state->lpdispersion >
	GET_DEFAULT(outlists->llk,1,0) + GET_DEFAULT(outlists->lpZ,1,0) + GET_DEFAULT(outlists->lpLV,1,0) + 
	GET_DEFAULT(outlists->lpcoef,1,0) + GET_DEFAULT(outlists->lpRE,1,0) + GET_DEFAULT(outlists->lpREV,1,0) +
	GET_DEFAULT(outlists->lpdispersion,1,0) )
      ERGMM_MCMC_store_iteration(1,model,cur->state,setting,outlists);

    /* every interval save the results */
    if((iter % setting->interval) == 0){
      pos = (iter/setting->interval-1)+ERGMM_OUTLISTS_RESERVE;

      // Store current iteration.
      ERGMM_MCMC_store_iteration(pos, model, cur->state, setting, outlists);

      // Acceptance rates.
      if(outlists->coef_rate){
        outlists->coef_rate[pos] = (double) ((double)n_accept_b)/((double)setting->interval);
      }
      if(outlists->Z_rate_move){
	outlists->Z_rate_move[pos] = (double) ((double)n_accept_z)/((double)setting->interval*model->verts); 
      }

      n_accept_z=0; 
      n_accept_b=0; 
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
  if(outlists->lpdispersion)
    outlists->lpdispersion[pos] = par->lpdispersion;

  // Covariate coefficients.
  Rpack_dvectors(par->coef,model->coef,outlists->coef+pos,setting->sample_size+ERGMM_OUTLISTS_RESERVE);


  if(model->latent){
    // Latent positions.
    Rpack_dmatrixs(par->Z,model->verts,model->latent,outlists->Z+pos,setting->sample_size+ERGMM_OUTLISTS_RESERVE);

      
    // Cluster-related.
    if(model->clusters>0){      
      // Cluster assignments.
      Rpack_ivectors((int *)par->Z_K,model->verts,outlists->Z_K+pos,setting->sample_size+ERGMM_OUTLISTS_RESERVE);
	  
      // Cluster means.
      Rpack_dmatrixs(par->Z_mean,model->clusters,model->latent,outlists->Z_mean+pos,setting->sample_size+ERGMM_OUTLISTS_RESERVE);
	  
      // Intracluster variances.
      Rpack_dvectors(par->Z_var,model->clusters,outlists->Z_var+pos,setting->sample_size+ERGMM_OUTLISTS_RESERVE);

      // Cluster probabilities.
      Rpack_dvectors(par->Z_pK,model->clusters,outlists->Z_pK+pos,setting->sample_size+ERGMM_OUTLISTS_RESERVE);

    }
    else
      outlists->Z_var[pos]=par->Z_var[0];
  }

  // Sender effects.
  if(par->sender){
    Rpack_dvectors(par->sender,model->verts,outlists->sender+pos,setting->sample_size+ERGMM_OUTLISTS_RESERVE);
    outlists->sender_var[pos] = par->sender_var;
  }

  // Receiver effects.
  if(par->receiver && !model->sociality){
    Rpack_dvectors(par->receiver,model->verts,outlists->receiver+pos,setting->sample_size+ERGMM_OUTLISTS_RESERVE);
    outlists->receiver_var[pos] = par->receiver_var;
  }      

  // Dispersion.
  if(model->dispersion){
    outlists->dispersion[pos] = par->dispersion;
  }      
      
}

