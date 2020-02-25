/*  File src/ergmm_probs.c in package latentnet, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2020 Statnet Commons
 */
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
#include "ergmm_latent_effects.h"
#include "ergmm_probs.h"

//#define VERBOSE 1
//#define SUPER_VERBOSE 1
//#define ALWAYS_RECOMPUTE_LLK 1

/* A wrapper for lp_Y_recompute to be called directly from R. */
void ERGMM_lp_Y_wrapper(int *n, int *p, int *d, int *latent_eff, int *family, int *res,
			int *dir, int *viY, double *vdY,
			int *iconsts, double *dconsts,
			double *vX, double *vZ,
			double *coef,
			double *sender, double *receiver,
			double *dispersion,
			int *vobserved_ties,
			double *llk){

  // This was added because newer versions of R no longer pass a 0-length vector as NULL, so we have to do it here.
  if(*p==0){
    vX=coef=NULL;
  }
  if(*d==0){
    vZ=NULL;
  }
  if(res[0]==0&&res[2]==0){
    sender=NULL;
  }
  if(res[1]==0){
    receiver=NULL;
  }

  unsigned int i,j,k;
  double **Z = vZ ? Runpack_dmatrix(vZ,*n,*d, NULL) : NULL;
  int **iY = viY ? Runpack_imatrix(viY, *n, *n, NULL) : NULL;
  double **dY = vdY ? Runpack_dmatrix(vdY, *n, *n, NULL) : NULL;
  unsigned int **observed_ties = (unsigned int **) (*vobserved_ties>=0 ? Runpack_imatrix(vobserved_ties,*n,*n,NULL) : NULL);
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
			    0,
			    res[2],
			    dispersion!=0?1:0,
			    *latent_eff ? ERGMM_MCMC_latent_eff[*latent_eff-1] : NULL
  };
  
  // Precompute the normalizing constant.
  ERGMM_MCMC_set_lp_Yconst[*family-1](&model);

  ERGMM_MCMC_Par params = {Z, // Z
			   coef, // coef
			   NULL, // Z_mean
			   NULL, // Z_var
			   NULL, // Z_pK			  
			   sender,
			   0, // sender_var
			   receiver,
			   0, // receiver_var
			   *dispersion,
			   NULL, // Z_K
			   0, // llk
			   NULL, // lpedge
			   0, // lpZ		  
			   0, // lpLV
			   0, // lpcoef
			   0, // lpRE
			   0  // lpdispersion
  };

  *llk = ERGMM_MCMC_lp_Y(&model,&params,FALSE);

  P_free_all();

  /* Memory freed by GC. */

  return;

}

/*R_INLINE*/ double ERGMM_MCMC_etaij(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par,unsigned int i,unsigned int j){
  double eta=0;
  unsigned int k;
  if(model->latent) eta+=model->latent_eff(par->Z[i],par->Z[j],model->latent);
  
  for(k=0;k<model->coef;k++) eta+=par->coef[k]*model->X[k][i][j];
  
  if(par->sender) eta+=par->sender[i];
  if(par->receiver) eta+=par->receiver[j];
  
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
  
  /* Do NOT add on log P(Z|means,vars,clusters) or log
     P(sender,receiver|vars) here.  This is a _log-likelihood_, so
     it's log P(Y|everything else), not log P(Y,Z|everything else).*/

  if(!update_lpedge) return(llk);
  else return(par->llk=llk);
}



double ERGMM_MCMC_lp_Y_diff(ERGMM_MCMC_Model *model, ERGMM_MCMC_MCMCState *cur){
  double llk_diff=0;
  unsigned int i,j,prop_i=PROP_NONE,prop_j=PROP_NONE;

  ERGMM_MCMC_Par *new=cur->prop,*old=cur->state;
  if(cur->prop_coef!=PROP_NONE || cur->prop_Z==PROP_ALL || cur->prop_RE==PROP_ALL)
    return(ERGMM_MCMC_lp_Y(model,new,TRUE)-old->llk);
  if(cur->prop_Z!=PROP_NONE) prop_i=prop_j=cur->prop_Z;
  else if(cur->prop_RE!=PROP_NONE){
    if(old->sender) prop_i=cur->prop_RE;
    if(old->receiver || model->sociality) prop_j=cur->prop_RE;
  }
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
  /* Do NOT add on log P(Z|means,vars,clusters) or log P(sender,receiver|vars) here.
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

  if(cur->prop_Z==PROP_ALL || cur->prop_LV!=PROP_NONE){
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


/* logp_RE gives log P(sender,receiver|sender_var,receiver_var)
 */
double ERGMM_MCMC_logp_RE(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par){
  unsigned int i;
  par->lpRE=0;
  for(i=0;i<model->verts;i++){
    if(par->sender) par->lpRE += dnorm(par->sender[i],0,sqrt(par->sender_var),1);
  }
  // Don't add receiver probability if the sender and receiver are "locked".
  if(par->receiver && !model->sociality)
    for(i=0;i<model->verts;i++)
      par->lpRE += dnorm(par->receiver[i],0,sqrt(par->receiver_var),1);

  return(par->lpRE);
}

double ERGMM_MCMC_logp_RE_diff(ERGMM_MCMC_Model *model, ERGMM_MCMC_MCMCState *cur){
  unsigned int i;
  double lpRE_diff=0;
  ERGMM_MCMC_Par *new=cur->prop,*old=cur->state;

  if(cur->prop_RE==PROP_ALL || cur->prop_REV!=PROP_NONE){
    return(ERGMM_MCMC_logp_RE(model,new)-old->lpRE);
  }
  else if(cur->prop_RE==PROP_NONE){
    new->lpRE=old->lpRE;
    return(0);
  }
  
  i=cur->prop_RE;

  if(new->sender) 
    lpRE_diff += (dnorm(new->sender[i],0,sqrt(new->sender_var),1)-
		  dnorm(old->sender[i],0,sqrt(old->sender_var),1));
  
  // Don't add receiver probability if the sender and receiver are "locked".
  if(new->receiver && !model->sociality)
    lpRE_diff += (dnorm(new->receiver[i],0,sqrt(new->receiver_var),1)-
		  dnorm(old->receiver[i],0,sqrt(old->receiver_var),1));

  new->lpRE=old->lpRE+lpRE_diff;
  return(lpRE_diff);
}

/* logp_latentvars gives P(Z_mean,Z_var|priors)
 * Note that it does NOT include probabilities involved in cluster assignments (might add later).
 */
double ERGMM_MCMC_logp_LV(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, ERGMM_MCMC_Priors *prior){
  unsigned int i,j;
  par->lpLV=0;
  if(model->clusters>0)
    for(i = 0; i < model->clusters; i++){
      for(j = 0; j < model->latent; j++)
	par->lpLV += dnorm(par->Z_mean[i][j],0,sqrt(prior->Z_mean_var),TRUE);
      par->lpLV+=dsclinvchisq(par->Z_var[i],prior->Z_var_df,prior->Z_var,TRUE);
    }
  else
    par->lpLV =dsclinvchisq(par->Z_var[0],prior->Z_var_df,prior->Z_var,TRUE);
  return(par->lpLV);
}

double ERGMM_MCMC_logp_LV_diff(ERGMM_MCMC_Model *model, ERGMM_MCMC_MCMCState *cur, ERGMM_MCMC_Priors *prior){
  if(cur->prop_LV==PROP_NONE){
    cur->prop->lpLV=cur->state->lpLV;
    return(0);
  }
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
  if(cur->prop_coef==PROP_NONE){
    cur->prop->lpcoef=cur->state->lpcoef;
    return(0);
  }
  return(ERGMM_MCMC_logp_coef(model,cur->prop,prior)-cur->state->lpcoef);
}

double ERGMM_MCMC_logp_REV(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, ERGMM_MCMC_Priors *prior){
  par->lpREV=0;
  if(par->sender) par->lpREV+=dsclinvchisq(par->sender_var,prior->sender_var_df,prior->sender_var,1);
  if(par->receiver && !model->sociality) par->lpREV+=dsclinvchisq(par->receiver_var,prior->receiver_var_df,prior->receiver_var,1);
  return(par->lpREV);
}

double ERGMM_MCMC_logp_REV_diff(ERGMM_MCMC_Model *model, ERGMM_MCMC_MCMCState *cur, ERGMM_MCMC_Priors *prior){
  if(cur->prop_REV==PROP_NONE){
    cur->prop->lpREV=cur->state->lpREV;
    return(0);
  }
  return(ERGMM_MCMC_logp_REV(model,cur->prop,prior)-cur->state->lpREV);
}

double ERGMM_MCMC_logp_dispersion(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, ERGMM_MCMC_Priors *prior){
  par->lpdispersion=0;
  if(model->dispersion) par->lpdispersion+=dsclinvchisq(par->dispersion,prior->dispersion_var_df,prior->dispersion_var,1);
  return(par->lpdispersion);
}

double ERGMM_MCMC_logp_dispersion_diff(ERGMM_MCMC_Model *model, ERGMM_MCMC_MCMCState *cur, ERGMM_MCMC_Priors *prior){
  if(cur->prop_dispersion==PROP_NONE){
    cur->prop->lpdispersion=cur->state->lpdispersion;
    return(0);
  }
  return(ERGMM_MCMC_logp_dispersion(model,cur->prop,prior)-cur->state->lpdispersion);
}
