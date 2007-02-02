/****************************************************************************/
/*  Original Author: Susan Shortreed, susanms@stat.washington.edu           */
/*  Updated by: Jeremy Tantrum, tantrum@stat.washington.edu                 */
/*  Purpose: support functions for parameter estimation model 2             */
/*           All of this code is for an R function which is incorporated    */
/*           into the R ERGM package.                                       */
/****************************************************************************/
/*####                 General Notes good to know                      #####*/
/*#  Z is the matrix of positions                                           */
/*#  k is the  dimension of the latent space                                */
/*#  g is the number of actors/nodes in the graph                           */
/*#                                                                         */
/****************************************************************************/

#include <R.h>
#include <Rmath.h>
#include <math.h> 
#include "ergmm_structs.h"
#include "matrix_utils.h"
#include "ergmm_utils.h"
#include "mvnorm.h"
#include "wishart.h"
#include "ergmm_probs.h"
#include "ergmm_updaters.h"


#define FALSE 0
#define TRUE !0

//#define VERBOSE 1
//#define SUPER_VERBOSE 1
//#define ALWAYS_RECOMPUTE_LLK 1


/*
  Initiates a Metropolis-Hasting step, by signaling what is about to be proposed.
  *** It MUST be called BEFORE the proposals are made. ***
  At the moment, it's not possible to propose LV without Z and REV without RE;
  this may change in the future. Gibbs-updating is OK, as long as the affected
  probabilities are updated.
 */
void ERGMM_MCMC_propose(ERGMM_MCMC_Model *model, ERGMM_MCMC_MCMCState *cur, unsigned int Z, unsigned int RE, unsigned int coef, unsigned int LV, unsigned int REV){
  // If the state has been Gibbs-updated, it means prop is inconsistent.
  if(cur->after_Gibbs) copy_MCMC_Par(model,cur->state,cur->prop);
  cur->after_Gibbs=FALSE;

  /* If we want to make a proposal for more than one vertex's Z, we have to
     Recalculate everything from scratch.
   */
  if(Z!=PROP_NONE && cur->state->Z){
    if(cur->prop_Z!=PROP_NONE && cur->prop_Z!=Z) cur->prop_Z=PROP_ALL;
    else cur->prop_Z=Z;
  }
  
  // Ditto RE...
  if(RE!=PROP_NONE && (cur->state->sender || cur->state->receiver )){
    if(cur->prop_RE!=PROP_NONE && cur->prop_RE!=RE) cur->prop_RE=PROP_ALL;
    else cur->prop_RE=RE;
  }

  if(coef!=PROP_NONE) cur->prop_coef=PROP_ALL;

  if(LV!=PROP_NONE && cur->state->Z){
    cur->prop_LV=PROP_ALL;
    cur->prop_Z=PROP_ALL;
  }
  if(REV!=PROP_NONE && (cur->state->sender || cur->state->receiver)){
    cur->prop_REV=PROP_ALL;
    cur->prop_RE=PROP_ALL;
  }

  /* Also, if we want to change latent position of one vertex but a random
     effect of another. */
  if(cur->prop_RE!=PROP_NONE && cur->prop_Z!=PROP_NONE && cur->prop_RE != cur->prop_Z)
    cur->prop_RE=cur->prop_Z = PROP_ALL;
}



void ERGMM_MCMC_prop_end(ERGMM_MCMC_Model *model, ERGMM_MCMC_MCMCState *cur,
	      ERGMM_MCMC_Par *new, ERGMM_MCMC_Par *old, unsigned int copy_lpedge){
  unsigned int i,j;

  // Overwrite old parameter values with new ones.
  switch(cur->prop_Z){
  case PROP_ALL:
    copy_dmatrix(new->Z,old->Z,model->verts,model->latent); break;
  case PROP_NONE: break;
  default:
    copy_dvector(new->Z[cur->prop_Z],old->Z[cur->prop_Z],model->latent); break;
  }

  switch(cur->prop_RE){
  case PROP_ALL:
    if(new->sender)
      copy_dvector(new->sender,old->sender,model->verts); 
    if(new->receiver && !model->sociality)
      copy_dvector(new->receiver,old->receiver,model->verts);
    break;
  case PROP_NONE: break;
  default:
    if(new->sender)
      old->sender[cur->prop_RE]=new->sender[cur->prop_RE];
    if(new->receiver && !model->sociality)
      old->receiver[cur->prop_RE]=new->receiver[cur->prop_RE];
  }
  
  if(cur->prop_coef==PROP_ALL)
    copy_dvector(new->coef,old->coef,model->coef);

  if(cur->prop_LV==PROP_ALL){
    if(new->Z_mu) copy_dmatrix(new->Z_mu,old->Z_mu,model->clusters,model->latent);
    if(new->Z_var) copy_dvector(new->Z_var,old->Z_var,model->clusters?model->clusters:1);
  }

  if(cur->prop_REV==PROP_ALL){
    if(new->sender) old->sender_var=new->sender_var;

    if(new->sender && model->sociality) old->receiver_var=new->sender_var;
    else if(new->receiver) old->receiver_var=new->receiver_var;      
  }

  // Update the lpedge matrix.
  if(copy_lpedge){
    if(cur->prop_Z==PROP_ALL || cur->prop_RE==PROP_ALL || cur->prop_coef==PROP_ALL){
      copy_dmatrix(new->lpedge,old->lpedge,model->verts,model->verts);
    }
    else if(cur->prop_Z!=PROP_NONE){
      i=cur->prop_Z;
      if(model->dir){
	copy_dvector(new->lpedge[i],old->lpedge[i],model->verts);
	for(j=0;j<i;j++) old->lpedge[j][i]=new->lpedge[j][i];
	for(j=i+1;j<model->verts;j++) old->lpedge[j][i]=new->lpedge[j][i];
      }
      else{
	copy_dvector(new->lpedge[i],old->lpedge[i],i);
	for(j=i+1;j<model->verts;j++) old->lpedge[j][i]=new->lpedge[j][i];
      }
    }
    else if(cur->prop_RE!=PROP_NONE){
      if(new->sender){
	i=cur->prop_RE;
	if(model->dir){
	  copy_dvector(new->lpedge[i],old->lpedge[i],model->verts);
	}
	else{
	  copy_dvector(new->lpedge[i],old->lpedge[i],i);
	  for(j=i+1;j<model->verts;j++) old->lpedge[j][i]=new->lpedge[j][i];
	}
      }
      /* Note that whether we update the "column" is determined by directedness
	 of the graph and presence of receiver effect, not whether it's locked
	 to sender effect. */
      if(new->receiver && model->dir){
	j=cur->prop_RE;
	for(i=0;i<model->verts;i++) old->lpedge[i][j]=new->lpedge[i][j];
      }
    }
  }

  // Finally, update probabilities and signal the end of proposal.
  old->llk=new->llk;

  if(cur->prop_Z!=PROP_NONE){
    old->lpZ=new->lpZ;
    cur->prop_Z=PROP_NONE;
  }

  if(cur->prop_RE!=PROP_NONE){
    old->lpRE=new->lpRE;
    cur->prop_RE=PROP_NONE;
  }

  if(cur->prop_coef!=PROP_NONE){
    old->lpcoef=new->lpcoef;
    cur->prop_coef=PROP_NONE;
  }

  if(cur->prop_LV!=PROP_NONE){
    old->lpLV=new->lpLV;
    cur->prop_LV=PROP_NONE;
  }

  if(cur->prop_REV!=PROP_NONE){
    old->lpREV=new->lpREV;
    cur->prop_REV=PROP_NONE;
  }
}

/* update Z and RE one vertex at a time */
unsigned int ERGMM_MCMC_Z_RE_up(ERGMM_MCMC_Model *model, ERGMM_MCMC_Priors *prior, ERGMM_MCMC_MCMCState *cur,
	    ERGMM_MCMC_MCMCSettings *setting)
{
  double lr;
  unsigned int iord, i, j, change=0;
  ERGMM_MCMC_Par *par=cur->prop;

  // Generate the order in which the vertices will be updated.
  runifperm(model->verts,cur->update_order);

  for(iord=0;iord<model->verts;iord++){
    i=cur->update_order[iord];
    // Make the proposal...
    
    ERGMM_MCMC_propose(model,cur,i,i,PROP_NONE,PROP_NONE,PROP_NONE);
    if(model->latent){
      for(j=0;j<model->latent;j++){ 
	par->Z[i][j] = cur->state->Z[i][j] + rnorm(0,setting->Z_delta);
      }
    }
    if(par->sender){
      par->sender[i] += rnorm(0,setting->RE_delta);
    }

    if(par->receiver && !model->sociality){
      par->receiver[i] += rnorm(0,setting->RE_delta);
    }


    lr = ERGMM_MCMC_lp_Y_diff(model,cur)+ERGMM_MCMC_logp_Z_diff(model,cur)+ERGMM_MCMC_logp_RE_diff(model,cur);

    if( runif(0.0,1.0) < exp(lr) ){
      change++;
      ERGMM_MCMC_accept(model,cur);
    }
    else{
      ERGMM_MCMC_reject(model,cur);
    }
  }

  return(change);
}

/* updates coef, scale of Z, and shifts the random effects; also translates Z */
unsigned int ERGMM_MCMC_coef_up_scl_tr_Z_shift_RE(ERGMM_MCMC_Model *model,  ERGMM_MCMC_Priors *prior, ERGMM_MCMC_MCMCState *cur,
			     ERGMM_MCMC_MCMCSettings *setting){  
  unsigned int i;
  double lr, h, dens_change=0;
  
  ERGMM_MCMC_Par *par=cur->prop;
  
  // Signal the proposal (of everything).
  ERGMM_MCMC_propose(model,cur,PROP_ALL,PROP_ALL,PROP_ALL,PROP_ALL,PROP_NONE);

  if(model->latent){  
    // Propose to scale Z.
    // Note that log P(mu,sigma) is changed.
    
    h = rlnorm(0,setting->Z_scl_delta);

    for(i=0;i<model->latent;i++)
      cur->tr_by[i]=rnorm(0,setting->Z_tr_delta);

    dmatrix_scale_by(par->Z,model->verts,model->latent,h);
    latentpos_translate(par->Z,model->verts,model->latent,cur->tr_by);
    if(model->clusters){
      dmatrix_scale_by(par->Z_mu,model->clusters,model->latent,h);
      latentpos_translate(par->Z_mu,model->clusters,model->latent,cur->tr_by);
      dvector_scale_by(par->Z_var,model->clusters,h*h);
    }
    else
      dvector_scale_by(par->Z_var,1,h*h);
  }

  // Propose coef.
  for(i=0;i< model->coef;i++){
    h=rnorm(0.0,setting->coef_delta[i]);
    par->coef[i] += h;
    dens_change += setting->X_means[i]*h;
  }

  // Propose to shift random effects.
  /* The proposal here effectively counteracts the change in density due to the proposal
     of coef.
     The SD of the change in density after all is said and done is settings->RE_shift_delta.
  */
  // It may make sense to also propose a jump in random effects variance here.
  
  //  dens_change=0;
  if(par->sender){
    if(model->sociality || !par->receiver){
      randeff_translate(par->sender,model->verts,
			(dens_change*ERGMM_MCMC_COEF_RE1_12
			 +rnorm(0.0,setting->RE_shift_delta)*ERGMM_MCMC_COEF_RE1_22));
    }else{
      h=rnorm(0.0,setting->RE_shift_delta);
      randeff_translate(par->sender,model->verts,
			(dens_change*ERGMM_MCMC_COEF_RE2_12
			 +h*ERGMM_MCMC_COEF_RE2_22));
    }
  }
  
  if(par->receiver && !model->sociality){
    if(par->sender)
      randeff_translate(par->receiver,model->verts,
			(dens_change*ERGMM_MCMC_COEF_RE2_13
			 +h*ERGMM_MCMC_COEF_RE2_23
			 +rnorm(0.0,setting->RE_shift_delta)*ERGMM_MCMC_COEF_RE2_33));
    else
      randeff_translate(par->receiver,model->verts,
			(dens_change*ERGMM_MCMC_COEF_RE1_12
			 +rnorm(0.0,setting->RE_shift_delta)*ERGMM_MCMC_COEF_RE1_22));
  }

  /* Calculate the log-likelihood-ratio.
     Note that even functions that don't make sense in context
     (e.g. logp_Z for a non-latent-space model) are safe to call and return 0). */
  lr = (ERGMM_MCMC_lp_Y_diff(model,cur)
	+ERGMM_MCMC_logp_coef_diff(model,cur,prior)
	+ERGMM_MCMC_logp_Z_diff(model,cur)
	+ERGMM_MCMC_logp_LV_diff(model,cur,prior)
	+ERGMM_MCMC_logp_RE_diff(model,cur));
  
  
  /* accept or reject */
  if( runif(0.0,1.0) < exp(lr) ){
    ERGMM_MCMC_accept(model,cur);
    return(1);
  }
  else{
    ERGMM_MCMC_reject(model,cur);
    return(0);
  }
}

void ERGMM_MCMC_CV_up(ERGMM_MCMC_Model *model, ERGMM_MCMC_Priors *prior, ERGMM_MCMC_MCMCState *cur){
  double S_hat, useSig,sum,temp;
  unsigned int i,j,loop1;

  ERGMM_MCMC_Par *par=cur->state;

  // Signal that this is a Gibbs-sampling update, so prop will need be caught up.
  cur->after_Gibbs=TRUE;

  // Reassign clusters.
  for(i=0;i<model->verts;i++){
    sum = 0;
    for(j=0;j<model->clusters;j++){
      temp=dindnormmu(model->latent,par->Z[i],par->Z_mu[j],sqrt(par->Z_var[j]),FALSE);
      if(j>0)
	cur->pK[j] = cur->pK[j-1] + par->Z_pK[j] * temp;
      else 
	cur->pK[j] = par->Z_pK[j] * temp;
    }
    
    temp = runif(0.0,1.0);
    j = 0;
    while(cur->pK[j]/cur->pK[model->clusters-1] < temp)
      j++;
    par->Z_K[i] = j + 1;
  }

  // Count the number of vertices in each cluster.
  for(i=0;i<model->clusters;i++)
    cur->n[i] = 0;
  for(i=0;i<model->verts;i++)
    cur->n[par->Z_K[i] - 1]++;
  
  //rdirichlet for epsilon
  for(i=0;i<model->clusters;i++){
    par->Z_pK[i] = (double)(cur->n[i]+prior->clust_dirichlet);
  }
  rdirich(model->clusters,par->Z_pK);
  
  // Compute cluster means
  init_dmatrix(cur->Z_bar,model->clusters,model->latent,0.0);
  for(i=0;i<model->verts;i++)
    for(j=0;j<model->latent;j++)
      cur->Z_bar[par->Z_K[i]-1][j] += par->Z[i][j]/cur->n[par->Z_K[i]-1];

  // invchisq for sigma
  for(i=0;i<model->clusters;i++){
    S_hat = 0.0;
    for(j=0;j<model->verts;j++)
      if((par->Z_K[j] - 1) == i)
	for(loop1=0;loop1<model->latent;loop1++)
	  S_hat += (par->Z[j][loop1] - par->Z_mu[i][loop1]) * (par->Z[j][loop1] - par->Z_mu[i][loop1]);

    par->Z_var[i] = rsclinvchisq(cur->n[i]*model->latent + prior->Z_var_df,
				 (prior->Z_var*prior->Z_var_df + S_hat)/
				 (cur->n[i]*model->latent + prior->Z_var_df));
    
  }

  //mvrnorm for mumat
  for(i=0;i<model->clusters;i++){
    for(j=0;j<model->latent;j++)
      cur->Z_bar[i][j] = cur->n[i] * cur->Z_bar[i][j]/(cur->n[i]+ par->Z_var[i]/prior->Z_mu_var);
    useSig = par->Z_var[i]/(cur->n[i] + par->Z_var[i]/prior->Z_mu_var);
    for(j=0;j<model->latent;j++)
      par->Z_mu[i][j] = rnorm(cur->Z_bar[i][j],sqrt(useSig));
  }
  // The following functions update par->lpZ and par->lpLV; they are NOT frivolous.
  ERGMM_MCMC_logp_Z(model,par);
  ERGMM_MCMC_logp_LV(model,par,prior);
}

/* Update latent variance.
 */
void ERGMM_MCMC_LV_up(ERGMM_MCMC_Model *model, ERGMM_MCMC_Priors *prior, ERGMM_MCMC_MCMCState *cur){
  unsigned int j,loop1;
  ERGMM_MCMC_Par *par=cur->state;

  // Signal that this is a Gibbs-sampling update, so prop will need be caught up.
  cur->after_Gibbs=TRUE;

  // invchisq for sigma
  double S_hat=0;
  for(j=0;j<model->verts;j++)
    for(loop1=0;loop1<model->latent;loop1++)
      S_hat += par->Z[j][loop1] * par->Z[j][loop1];

  par->Z_var[0] = rsclinvchisq(model->verts*model->latent + prior->Z_var_df,
			       (prior->Z_var*prior->Z_var_df + S_hat)/
			       (model->verts*model->latent + prior->Z_var_df));

  // The following functions update par->lpZ and par->lpLV; they are NOT frivolous.
  ERGMM_MCMC_logp_Z(model,par);
  ERGMM_MCMC_logp_LV(model,par,prior);
}

void ERGMM_MCMC_REV_up(ERGMM_MCMC_Model *model, ERGMM_MCMC_Priors *prior, ERGMM_MCMC_MCMCState *cur){
  double S_hat;
  unsigned int i;
  ERGMM_MCMC_Par *par=cur->state;

  // Signal that this is a Gibbs-sampling update, so prop will need be caught up.
  cur->after_Gibbs=TRUE;

  //Inverse-chisq (Wishart in the future?) for sender and receiver effect variances.
  if(par->sender){
    S_hat = 0.0;
    for(i=0;i<model->verts;i++)
      S_hat += par->sender[i] * par->sender[i];
    par->sender_var = rsclinvchisq(model->verts + prior->sender_var_df,
				   (prior->sender_var*prior->sender_var_df +  S_hat)/
				   (model->verts + prior->sender_var_df));
  }


  if(par->receiver && !model->sociality){
    S_hat = 0.0;
    for(i=0;i<model->verts;i++)
      S_hat += par->receiver[i] * par->receiver[i];
    par->receiver_var = rsclinvchisq(model->verts + prior->receiver_var_df,
				     (prior->receiver_var*prior->receiver_var_df + S_hat)/
				     (model->verts + prior->receiver_var_df));
  }
  else par->receiver_var=par->sender_var;

  // The following function updates par->lpRE; it is NOT frivolous.
  ERGMM_MCMC_logp_RE(model,par);
  ERGMM_MCMC_logp_REV(model,par,prior);
}
