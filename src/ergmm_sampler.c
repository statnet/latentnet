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
			int *viY,
			double *vdY,
			int *family,
			int *iconsts,
			double *dconsts,

			double *vX,
			  
			double *llk_mcmc,
			double *lpZ_mcmc,
			double *lpcoef_mcmc,
			double *lpLV_mcmc,
			   
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
			double *Z_rate_move_all,

			int *Z_K_mcmc,
			double *Z_pK_mcmc,
			double *Z_mean_mcmc,
			double *Z_var_mcmc,
			  
			double *coef_start,
			double *coef_prior_mean,
			double *coef_var,
			double *coef_mcmc,
			double *coef_rate, 
			  
			int *vobserved_ties,
			double *deltas){
  int i=0,j=0,k=0;
  double **Z_start = vZ_start ? Runpack_dmatrix(vZ_start,*n,*d, NULL) : NULL;
  double **Z_mean_start = vZ_mean_start ? Runpack_dmatrix(vZ_mean_start,*G,*d,NULL) : NULL;
  int **iY = viY ? Runpack_imatrix(viY, *n, *n, NULL):NULL;
  double **dY = vdY ? Runpack_dmatrix(vdY, *n, *n, NULL):NULL;
  unsigned int **observed_ties = vobserved_ties ? (unsigned int **) Runpack_imatrix(vobserved_ties,*n,*n,NULL) : NULL;
  double ***X = d3array(*p,*n,*n);

  // set up all of the covariate matrices if covariates are involed 
  // if p=0 (ie no covariates then these next two loops will do nothing)
  //

  for(k=0;k<*p;k++){
    for(i=0;i<*n;i++){
      for(j=0;j<*n;j++){
	X[k][i][j] = vX[ k*(*n)*(*n) + i*(*n) + j ];
      }
    }
  }

  /* R function enabling uniform RNG */
  GetRNGstate();
 

  ERGMM_MCMC_init(*samples_stored, *interval,

		  *n,*p,
		  d ? *d : 0,
		  *G,
		 
		  *dir,iY,dY,
		  family ? *family-1 : 0,iconsts,dconsts,
		  X,

		  llk_mcmc, lpZ_mcmc, lpcoef_mcmc, lpLV_mcmc,
		    
		  Z_start, 
		  Z_pK_start,Z_mean_start,Z_var_start,(unsigned int *)Z_K_start,
		  Z_var_prior? *Z_var_prior : 0,
		  Z_mean_prior_var ? *Z_mean_prior_var : 0,
		  Z_pK_prior ? *Z_pK_prior : 0,
		  Z_var_prior_df ? *Z_var_prior_df : 0,
		  Z_mcmc, Z_rate_move, Z_rate_move_all, Z_K_mcmc, Z_pK_mcmc, Z_mean_mcmc, Z_var_mcmc,

		  coef_start,
		  coef_mcmc, coef_rate,    
		  coef_prior_mean, coef_var,

		  observed_ties,
		  deltas[0],deltas[1],deltas[2],
		  deltas+COEF_DELTA_START);

  PutRNGstate();
  P_free_all();
  return;
}

/* Initializes the MCMC sampler and allocates memory. */
void ERGMM_MCMC_init(unsigned int samples_stored, unsigned int interval, 

		     unsigned int n,
		     unsigned int p, unsigned int d, unsigned int G,

		     unsigned int dir, int **iY, double **dY,

		     unsigned int family, int *iconsts, double *dconsts,

		     double ***X,

		     double *llk_mcmc, double *lpZ_mcmc, double *lpcoef_mcmc, double *lpLV_mcmc,

		     double **Z_start,
		     double *Z_pK_start, double **Z_mean_start, double *Z_var_start, unsigned int *Z_K_start,
		     double Z_var_prior, double Z_mean_prior_var, double Z_pK_prior,
		     double Z_var_prior_df,
		     double *Z_mcmc, double *Z_rate_move, double *Z_rate_move_all, int *K_mcmc,
		     double *Z_pK_mcmc,
		     double *Z_mean_mcmc, double *Z_var_mcmc,

		     double *coef_start,
		     double *coef_mcmc, double *coef_rate, 
		     double *coef_prior_mean, double *coef_var, 

		     unsigned int **observed_ties,

		     double Z_delta, double Z_tr_delta, double Z_scl_delta,
		     double *coef_delta)
{
  unsigned int i,j,k,n_observed;
  double X_sum;

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
			    G // clusters
			    };
  ERGMM_MCMC_set_lp_Yconst[family](&model);

  ERGMM_MCMC_MCMCSettings setting = {Z_delta,Z_tr_delta,Z_scl_delta,
				     coef_delta,
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

  ERGMM_MCMC_Priors prior = {Z_mean_prior_var, // Z_mean_var
			     Z_var_prior, // Z_var
			     Z_var_prior_df, // a.k.a. Z_var_df (I hope)
			     coef_prior_mean,
			     coef_var,
			     Z_pK_prior};
  
  ERGMM_MCMC_Par state = {Z_start, // Z
			  coef_start, // coef
			  Z_mean_start, // Z_mean
			  Z_var_start, // Z_var
			  Z_pK_start, // Z_pK			  
			  Z_K_start, // Z_K
			  0, // llk
			  dmatrix(model.verts,model.verts), // lpedge
			  0, // lpZ		  
			  0, // lpLV
			  0, // lpcoef
  };

  ERGMM_MCMC_Par prop = {model.latent ? dmatrix(model.verts,model.latent):NULL, // Z
			 model.coef ? dvector(model.coef):NULL, // coef
			 model.clusters ? dmatrix(model.clusters,model.latent):NULL, // Z_mean
			 model.latent ? dvector(model.clusters?model.clusters:1):NULL, // Z_var
			 model.clusters ? dvector(model.clusters):NULL, // Z_pK
			 state.Z_K, // prop.Z_K === state.Z_K
			 0, // llk
			 dmatrix(model.verts,model.verts), // lpedge
			 0, // lpZ
			 0, // lpLV
			 0, // lpcoef
  };

  ERGMM_MCMC_MCMCState start = {&state,
				&prop,
				model.clusters ? dmatrix(model.clusters,model.latent) : NULL, // Z_bar
				model.latent ? dvector(model.latent) : NULL, // tr_by
				model.clusters ? dvector(model.clusters): NULL, // pK
				model.clusters ? (unsigned int *) ivector(model.clusters) : NULL, // n
				PROP_NONE, // prop_Z
				PROP_NONE, // prop_coef
				PROP_NONE, // prop_LV
				FALSE, // after_Gibbs
				model.latent ? (unsigned int *) ivector(model.verts) : NULL // update_order
  };
  
  ERGMM_MCMC_ROutput outlists = {llk_mcmc, lpZ_mcmc, lpcoef_mcmc, lpLV_mcmc,
				 Z_mcmc, Z_rate_move, Z_rate_move_all,
				 coef_mcmc,coef_rate,
				 Z_mean_mcmc,Z_var_mcmc,Z_pK_mcmc,
				 K_mcmc};

  if(model.clusters>0)
    for(i=0;i<model.verts;i++)
      start.n[state.Z_K[i] - 1]++;

  // Initialize the log-probabilities.
  state.llk = ERGMM_MCMC_lp_Y(&model, &state, TRUE);
  if(model.latent) ERGMM_MCMC_logp_Z(&model, &state);
  if(state.coef) ERGMM_MCMC_logp_coef(&model, &state, &prior);
  if(model.latent) ERGMM_MCMC_logp_LV(&model, &state, &prior);
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

    n_accept_z += ERGMM_MCMC_Z_up(model, prior, cur, setting); 

    if(model->latent){

      //n_accept_transl_z += translate_Z(model,prior,cur,setting);

      // Update cluster parameters (they are separated from data by Z, so plain Gibbs).
      // Note that they are also updated in coef_up_scl_tr_Z.
      if(model->clusters>0)
	ERGMM_MCMC_CV_up(model,prior,cur);
      else
	ERGMM_MCMC_LV_up(model,prior,cur);
    }

    /* Update coef given this new value of Z and conditioned on everything else.
       Also propose to scale Z.
    */
    
    if( ERGMM_MCMC_coef_up_scl_tr_Z(model,prior,cur,setting) ){
      n_accept_b++;
    }  

    // If we have a new MLE (actually, the highest likelihood encountered to date), store it.
    if( cur->state->llk > GET_DEFAULT(outlists->llk,0,0) ) ERGMM_MCMC_store_iteration(0,model,cur->state,setting,outlists);

    // If we have a new posterior mode (actually, the highest joint density of all variables but K observed to date), store it.
    if( cur->state->llk + cur->state->lpZ + cur->state->lpLV + 
	cur->state->lpcoef >
	GET_DEFAULT(outlists->llk,1,0) + GET_DEFAULT(outlists->lpZ,1,0) + GET_DEFAULT(outlists->lpLV,1,0) + 
	GET_DEFAULT(outlists->lpcoef,1,0) )
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
  if(outlists->lpLV)
    outlists->lpLV[pos] = par->lpLV;

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
}

