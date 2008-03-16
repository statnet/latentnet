/***********************************************************************/
/* Functions to update individual ERGMM parameters from the posterior. */
/***********************************************************************/

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

/*
  Initiates a Metropolis-Hasting step, by signaling what is about to be proposed.
  *** It MUST be called BEFORE the proposals are made. ***
  Gibbs-updating is OK, as long as the affected probabilities are updated.
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
  }
  if(REV!=PROP_NONE && (cur->state->sender || cur->state->receiver)){
    cur->prop_REV=PROP_ALL;
  }

  /* Also, if we want to change latent position of one vertex but a random
     effect of another. */
  if(cur->prop_RE!=PROP_NONE && cur->prop_Z!=PROP_NONE && cur->prop_RE != cur->prop_Z)
    cur->prop_RE=cur->prop_Z = PROP_ALL;
}

/*
  Finalizes and cleans up after an MH proposal, copying and overwriting as needed.
  Usually not called directly, but through macros ERGMM_MCMC_accept and 
  ERGMM_MCMC_reject.
*/

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
    if(new->Z_mean) copy_dmatrix(new->Z_mean,old->Z_mean,model->clusters,model->latent);
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

  if(cur->prop_Z!=PROP_NONE){
    old->llk=new->llk;
    old->lpZ=new->lpZ;
    cur->prop_Z=PROP_NONE;
  }

  if(cur->prop_RE!=PROP_NONE){
    old->llk=new->llk;
    old->lpRE=new->lpRE;
    cur->prop_RE=PROP_NONE;
  }

  if(cur->prop_coef!=PROP_NONE){
    old->llk=new->llk;
    old->lpcoef=new->lpcoef;
    cur->prop_coef=PROP_NONE;
  }

  if(cur->prop_LV!=PROP_NONE){
    old->lpLV=new->lpLV;
    // A jump in latent space variables can also affect the probability of Z.
    old->lpZ=new->lpZ;
    cur->prop_LV=PROP_NONE;
  }

  if(cur->prop_REV!=PROP_NONE){
    old->lpREV=new->lpREV;
    // Ditto RE...
    old->lpRE=new->lpRE;
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

    lr = (ERGMM_MCMC_lp_Y_diff(model,cur)
	  +ERGMM_MCMC_logp_Z_diff(model,cur)
	  +ERGMM_MCMC_logp_RE_diff(model,cur)
	  );
	  
    if( setting->accept_all || runif(0.0,1.0) < exp(lr) ){
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
unsigned int ERGMM_MCMC_coef_up_scl_Z_shift_RE(ERGMM_MCMC_Model *model,  ERGMM_MCMC_Priors *prior, ERGMM_MCMC_MCMCState *cur,
					       ERGMM_MCMC_MCMCSettings *setting){
  double acc_adjust=0;
  ERGMM_MCMC_Par *par=cur->prop;
  
  // Signal the proposal (of everything).
  ERGMM_MCMC_propose(model,cur,PROP_ALL,PROP_ALL,PROP_ALL,PROP_ALL,PROP_NONE);

  for(unsigned int j=0; j<setting->group_prop_size; j++) 
    cur->deltas[j] = 0;

  for(unsigned int i=0; i<setting->group_prop_size; i++){
    double delta = rnorm(0,1);
    for(unsigned int j=0; j<setting->group_prop_size; j++){
      cur->deltas[j] += setting->group_deltas[i][j]*delta;
    }
  }
  
  unsigned int prop_pos=0;
  
  // Propose coef.
  for(unsigned int i=0;i< model->coef;i++){
    par->coef[i] += cur->deltas[prop_pos++];
  }

  if(model->latent){  
    // Propose to scale Z.
    // Note that log P(mu,sigma) is changed.

    // Grab the scaling factor.
    double logh=cur->deltas[prop_pos++];

    // Every scaling proposal has an accompanying adjustment for the log-probability
    // ratios.
    dmatrix_scale_by(par->Z,model->verts,model->latent,exp(logh));
    acc_adjust+=model->verts*model->latent*logh;

    if(model->clusters){
      dmatrix_scale_by(par->Z_mean,model->clusters,model->latent,exp(logh));
      acc_adjust+=model->clusters*model->latent*logh;

      dvector_scale_by(par->Z_var,model->clusters,exp(2*logh));
      acc_adjust+=model->clusters*2*logh;
    }else{
      dvector_scale_by(par->Z_var,1,exp(2*logh));
      acc_adjust+=2*logh;
    }
  }

  // Propose to shift random effects.
  if(par->sender){
    for(unsigned int i=0; i<model->verts; i++)
      par->sender[i]+=cur->deltas[prop_pos++];
  }
  
  if(par->receiver && !model->sociality){
    for(unsigned int i=0; i<model->verts; i++)
      par->receiver[i]+=cur->deltas[prop_pos++];
  }

  /* Calculate the log-likelihood-ratio.
     Note that even functions that don't make sense in context
     (e.g. logp_Z for a non-latent-space model) are safe to call and return 0). */

  double lr = (+ERGMM_MCMC_lp_Y_diff(model,cur)
	       +ERGMM_MCMC_logp_coef_diff(model,cur,prior)
	       +ERGMM_MCMC_logp_Z_diff(model,cur)
	       +ERGMM_MCMC_logp_LV_diff(model,cur,prior)
	       +ERGMM_MCMC_logp_RE_diff(model,cur)
	       +acc_adjust
	       );
  
  /* accept or reject */
  if( setting->accept_all || runif(0.0,1.0) < exp(lr) ){
    ERGMM_MCMC_accept(model,cur);
    return(1);
  }
  else{
    ERGMM_MCMC_reject(model,cur);
    return(0);
  }
}

/*
  Update "clustering variables" (mean, variance, and assignments).
*/

void ERGMM_MCMC_CV_up(ERGMM_MCMC_Model *model, ERGMM_MCMC_Priors *prior, ERGMM_MCMC_MCMCState *cur){
  double S_hat, useSig,temp;
  unsigned int i,j,loop1;

  ERGMM_MCMC_Par *par=cur->state;

  // Signal that this is a Gibbs-sampling update, so prop will need be caught up.
  cur->after_Gibbs=TRUE;

  // Reassign clusters.
  for(i=0;i<model->verts;i++){
    for(j=0;j<model->clusters;j++){
      temp=dindnormmu(model->latent,par->Z[i],par->Z_mean[j],sqrt(par->Z_var[j]),FALSE);
      if(j>0)
	cur->pK[j] = cur->pK[j-1] + par->Z_pK[j] * temp;
      else 
	cur->pK[j] = par->Z_pK[j] * temp;
    }
    
    temp = runif(0.0,1.0);
    for(j=0; cur->pK[j]/cur->pK[model->clusters-1] < temp; j++); // NOTE: this is not a bug; the for loop is there to find the right j!
    par->Z_K[i] = j + 1;
  }

  // Count the number of vertices in each cluster.
  for(i=0;i<model->clusters;i++)
    cur->n[i] = 0;
  for(i=0;i<model->verts;i++)
    cur->n[par->Z_K[i] - 1]++;
  
  // Dirichlet for epsilon
  for(i=0;i<model->clusters;i++){
    par->Z_pK[i] = (double)(cur->n[i]+prior->Z_pK);
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
	  S_hat += (par->Z[j][loop1] - par->Z_mean[i][loop1]) * (par->Z[j][loop1] - par->Z_mean[i][loop1]);

    par->Z_var[i] = rsclinvchisq(cur->n[i]*model->latent + prior->Z_var_df,
				 (prior->Z_var*prior->Z_var_df + S_hat)/
				 (cur->n[i]*model->latent + prior->Z_var_df));
    
  }

  // Z_mean
  for(i=0;i<model->clusters;i++){
    for(j=0;j<model->latent;j++)
      cur->Z_bar[i][j] = cur->n[i] * cur->Z_bar[i][j]/(cur->n[i]+ par->Z_var[i]/prior->Z_mean_var);
    useSig = par->Z_var[i]/(cur->n[i] + par->Z_var[i]/prior->Z_mean_var);
    for(j=0;j<model->latent;j++)
      par->Z_mean[i][j] = rnorm(cur->Z_bar[i][j],sqrt(useSig));
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
