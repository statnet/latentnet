/******************************************************/
/* Functions to compute probabilities and likelihoods */
/******************************************************/



#include <R.h>
#include <Rmath.h>
#include <math.h> 
#include "matrix_utils.h"
#include "ergmm_utils.h"
#include "mvnorm.h"
#include "wishart.h"
#include "ergmm_families.h"
#include "ergmm_probs.h"

#define FALSE 0
#define TRUE !0

//#define VERBOSE 1
//#define SUPER_VERBOSE 1
//#define ALWAYS_RECOMPUTE_LLK 1

/* A wrapper for lp_Y_recompute to be called directly from R. */
void ERGMM_lp_Y_wrapper(int *n, int *p, int *d,
			int *dir, int *viY, double *vdY,
			int *family, int *iconsts, double *dconsts,
			double *vX, double *vZ,
			double *coef,
			int *vobserved_ties,
			double *llk){
  unsigned int i,j,k;
  double **Z = vZ ? Runpack_dmatrix(vZ,*n,*d, NULL) : NULL;
  int **iY = viY ? Runpack_imatrix(viY, *n, *n, NULL) : NULL;
  double **dY = vdY ? Runpack_dmatrix(vdY, *n, *n, NULL) : NULL;
  unsigned int **observed_ties = (unsigned int **) (vobserved_ties ? Runpack_imatrix(vobserved_ties,*n,*n,NULL) : NULL);
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
  
  ERGMM_MCMC_Model model = {*dir,
			    iY, // iY
			    dY, // dY
			    X, // X
			    observed_ties,
			    ERGMM_MCMC_lp_edge[*family-1],
			    ERGMM_MCMC_E_edge[*family-1],
			    0,
			    iconsts,
			    dconsts,
			    *n, // verts
			    *d, // latent
			    *p, // coef
			    0};
  
  // Precompute the normalizing constant.
  ERGMM_MCMC_set_lp_Yconst[*family](&model);

  ERGMM_MCMC_Par params = {Z, // Z
			   coef, // coef
			   NULL, // Z_mean
			   NULL, // Z_var
			   NULL, // Z_pK			  
			   NULL, // Z_K
			   0, // llk
			   NULL, // lpedge
			   0, // lpZ		  
			   0, // lpLV
			   0, // lpcoef
  };

  *llk = ERGMM_MCMC_lp_Y(&model,&params,FALSE);

  P_free_all();

  /* Memory freed by GC. */

  return;

}

R_INLINE double ERGMM_MCMC_etaij(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par,unsigned int i,unsigned int j){
  double eta=0;
  unsigned int k;
  if(model->latent) eta-=dvector_dist(par->Z[i],par->Z[j],model->latent);
  
  for(k=0;k<model->coef;k++) eta+=par->coef[k]*model->X[k][i][j];
  
  return(eta);
}

/* recomputes log probability of the graph
   "own_lpedge" is used for debugging --- if set, lp_Y_recompute allocates its own
   lpedge matrix and deallocates it when done, leaving no traces of it having been run */
double ERGMM_MCMC_lp_Y(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int update_lpedge){
  double llk=model->lp_Yconst, lpedge;
  unsigned int i, j;

  /* actually calculating log-likelihood */
  // This is actually more readable than the nested ifs.
  /* I don't trust the compiler to optimize the simpler form properly.
     (It might not know to exchange if statements with loops.)
  */

  if(model->dir){
    for(i=0;i<model->verts;i++)
      for(j=0;j<model->verts;j++)
	if(IS_OBSERVABLE(model->observed_ties,i,j)){
	  llk+=lpedge=model->lp_edge(model,par,i,j);
	  if(update_lpedge) par->lpedge[i][j]=lpedge;
	}
  }
  else{
    for(i=0;i<model->verts;i++)
      for(j=0;j<=i;j++)
	if(IS_OBSERVABLE(model->observed_ties,i,j)){
	  llk+=lpedge=model->lp_edge(model,par,i,j);
	  if(update_lpedge) par->lpedge[i][j]=lpedge;
	}
  }
  
  /* Do NOT add on log P(Z|means,vars,clusters) here.  This is a _log-likelihood_, so
     it's log P(Y|everything else), not log P(Y,Z|everything else).*/

  if(!update_lpedge) return(llk);
  else return(par->llk=llk);
}



double ERGMM_MCMC_lp_Y_diff(ERGMM_MCMC_Model *model, ERGMM_MCMC_MCMCState *cur){
  double llk_diff=0;
  unsigned int i,j,prop_i=PROP_NONE,prop_j=PROP_NONE;

  ERGMM_MCMC_Par *new=cur->prop,*old=cur->state;
  if(cur->prop_coef!=PROP_NONE || cur->prop_Z==PROP_ALL)
    return(ERGMM_MCMC_lp_Y(model,new,TRUE)-old->llk);
  if(cur->prop_Z!=PROP_NONE) prop_i=prop_j=cur->prop_Z;
  else{
    new->llk=old->llk;
    return(0);
  }
  
  if(model->dir){ // Directed.
    if(prop_i!=PROP_NONE){
      i=prop_i;
      for(j=0; j<model->verts; j++){
	if(IS_OBSERVABLE(model->observed_ties,i,j)){
	  llk_diff+=new->lpedge[i][j]=model->lp_edge(model,new,i,j);
	  llk_diff-=old->lpedge[i][j];
	}
      }
    }
    if(prop_j!=PROP_NONE){
      j=prop_j;
      for(i=0; i<model->verts; i++){
	if(IS_OBSERVABLE(model->observed_ties,i,j)){
	  llk_diff+=new->lpedge[i][j]=model->lp_edge(model,new,i,j);
	  llk_diff-=old->lpedge[i][j];
	}
      }
    }
    // Don't doublecount edge (i,j).
    if(prop_i!=PROP_NONE && prop_j!=PROP_NONE){
      i=prop_i;
      j=prop_j;
      if(IS_OBSERVABLE(model->observed_ties,i,j)){
	llk_diff-=new->lpedge[i][j]=new->lpedge[i][j];
	llk_diff+=old->lpedge[i][j];
      }
    }
  }
  else{ // Undirected.
    if(prop_i!=PROP_NONE){
      i=prop_i;
      for(j=0; j<i; j++){
	if(IS_OBSERVABLE(model->observed_ties,i,j)){
	  llk_diff+=new->lpedge[i][j]=model->lp_edge(model,new,i,j);
	  llk_diff-=old->lpedge[i][j];
	}
      }
      for(j=i; j<model->verts; j++){
	if(IS_OBSERVABLE(model->observed_ties,j,i)){
	  llk_diff+=new->lpedge[j][i]=model->lp_edge(model,new,j,i);
	  llk_diff-=old->lpedge[j][i];
	}
      }
    }
  }

  new->llk=old->llk+llk_diff;
  /* Do NOT add on log P(Z|means,vars,clusters) here.
     This is a _log-likelihood_, so it's log P(Y|everything else), not log P(Y,Z|everything else).*/

  return(llk_diff);
}



/* logp_Z gives log P(Z|Z_mean,Z_K,Z_var) 
 */
double ERGMM_MCMC_logp_Z(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par){
  unsigned int i;
  par->lpZ=0;
  if(model->clusters>0)
    for(i=0;i<model->verts;i++)
      par->lpZ += dindnormmu(model->latent,par->Z[i],
			     par->Z_mean[par->Z_K[i]-1],sqrt(par->Z_var[par->Z_K[i]-1]),TRUE);
  else
    for(i=0;i<model->verts;i++)
      par->lpZ += diidnorm0(model->latent, par->Z[i], sqrt(par->Z_var[0]),TRUE);
  return(par->lpZ);
}

double ERGMM_MCMC_logp_Z_diff(ERGMM_MCMC_Model *model, ERGMM_MCMC_MCMCState *cur){
  unsigned int i;
  double lpZ_diff=0;
  ERGMM_MCMC_Par *new=cur->prop,*old=cur->state;

  if(cur->prop_Z==PROP_ALL){
    return(ERGMM_MCMC_logp_Z(model,new)-old->lpZ);
  }
  else if(cur->prop_Z==PROP_NONE){
    new->lpZ=old->lpZ;
    return(0);
  }
  
  i=cur->prop_Z;

  if(model->clusters>0)
    lpZ_diff=dindnormmu(model->latent,new->Z[i],new->Z_mean[new->Z_K[i]-1],sqrt(new->Z_var[new->Z_K[i]-1]),TRUE)-
      dindnormmu(model->latent,old->Z[i],old->Z_mean[old->Z_K[i]-1],sqrt(old->Z_var[old->Z_K[i]-1]),TRUE);
  else
    lpZ_diff=diidnorm0(model->latent, new->Z[i], sqrt(new->Z_var[0]),TRUE)-
      diidnorm0(model->latent, old->Z[i], sqrt(old->Z_var[0]),TRUE);

  new->lpZ=old->lpZ+lpZ_diff;
  return(lpZ_diff);
}


/* logp_LV gives P(Z_mean,Z_var|priors)
 * Note that it does NOT include probabilities involved in cluster assignments (might add later).
 */
double ERGMM_MCMC_logp_LV(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, ERGMM_MCMC_Priors *prior){
  unsigned int i,j;
  par->lpLV=0;
  if(model->clusters>0)
    for(i = 0; i < model->clusters; i++){
      for(j = 0; j < model->latent; j++)
	par->lpLV += dnorm(par->Z_mean[i][j],0,sqrt(prior->Z_mean_var),1);
      par->lpLV+=dsclinvchisq(par->Z_var[i],prior->Z_var_df,prior->Z_var,1);
    }
  else
    par->lpLV =dsclinvchisq(par->Z_var[0],prior->Z_var_df,prior->Z_var,1);
  return(par->lpLV);
}

double ERGMM_MCMC_logp_LV_diff(ERGMM_MCMC_Model *model, ERGMM_MCMC_MCMCState *cur, ERGMM_MCMC_Priors *prior){
  if(cur->prop_LV==PROP_NONE) return(0);
  return(ERGMM_MCMC_logp_LV(model,cur->prop,prior)-cur->state->lpLV);
}

/* logp_coef gives P(coef|priors)
 */
double ERGMM_MCMC_logp_coef(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, ERGMM_MCMC_Priors *prior){
  unsigned int i;
  par->lpcoef=0;
  for(i=0;i<model->coef;i++)
    par->lpcoef += dnorm(par->coef[i],prior->coef_mean[i],sqrt(prior->coef_var[i]),1);
  return(par->lpcoef);
}

double ERGMM_MCMC_logp_coef_diff(ERGMM_MCMC_Model *model, ERGMM_MCMC_MCMCState *cur, ERGMM_MCMC_Priors *prior){
  if(cur->prop_coef==PROP_NONE) return(0);
  return(ERGMM_MCMC_logp_coef(model,cur->prop,prior)-cur->state->lpcoef);
}
