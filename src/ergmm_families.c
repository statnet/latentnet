/*  File src/ergmm_families.c in package latentnet, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2018 Statnet Commons
 */
/****************************************************************/
/* Families of dyad weight distributions supported by latentnet */
/****************************************************************/

#include "ergmm_families.h"
#include "ergmm_probs.h"
#include <R.h>
#include <Rmath.h>
#include <math.h>


/* ERGMM_MCMC_lp_edge = unnormalized log-probability of a dyad. */
/* ERGMM_MCMC_set_lp_Yconst = (one-time) setting of the network's normalization constant. */
/* ERGMM_MCMC_E_edge = expected value of a dyad. */

/* Define "lookup tables" for families. */

const unsigned int ERGMM_MCMC_is_discrete[N_FAMILIES]={TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE};
const unsigned int ERGMM_MCMC_to_cont[N_FAMILIES]={3,4,5,3,4,5,6};

ERGMM_MCMC_lp_edge_t ERGMM_MCMC_lp_edge[N_FAMILIES]={
  ERGMM_MCMC_lp_edge_Bernoulli_logit,
  ERGMM_MCMC_lp_edge_binomial_logit,
  ERGMM_MCMC_lp_edge_Poisson_log,
  ERGMM_MCMC_lp_edge_Bernoulli_cont_logit,
  ERGMM_MCMC_lp_edge_binomial_cont_logit,
  ERGMM_MCMC_lp_edge_Poisson_cont_log,
  ERGMM_MCMC_lp_edge_normal_identity
};
  
ERGMM_MCMC_set_lp_Yconst_t ERGMM_MCMC_set_lp_Yconst[N_FAMILIES]={
  ERGMM_MCMC_set_lp_Yconst_Bernoulli_logit,
  ERGMM_MCMC_set_lp_Yconst_binomial_logit,
  ERGMM_MCMC_set_lp_Yconst_Poisson_log,
  ERGMM_MCMC_set_lp_Yconst_Bernoulli_cont_logit,
  ERGMM_MCMC_set_lp_Yconst_binomial_cont_logit,
  ERGMM_MCMC_set_lp_Yconst_Poisson_cont_log,
  ERGMM_MCMC_set_lp_Yconst_normal_identity
};
  
ERGMM_MCMC_E_edge_t ERGMM_MCMC_E_edge[N_FAMILIES]={
  ERGMM_MCMC_E_edge_Bernoulli_logit,
  ERGMM_MCMC_E_edge_binomial_logit,
  ERGMM_MCMC_E_edge_Poisson_log,
  ERGMM_MCMC_E_edge_Bernoulli_cont_logit,
  ERGMM_MCMC_E_edge_binomial_cont_logit,
  ERGMM_MCMC_E_edge_Poisson_cont_log,
  ERGMM_MCMC_E_edge_normal_identity
};


/*
  0 Bernoulli_logit
*/
double ERGMM_MCMC_lp_edge_Bernoulli_logit(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j){
  double eta=ERGMM_MCMC_etaij(model,par,i,j);
  return(model->iY[i][j]*eta-log1p(exp(eta)));
}

void ERGMM_MCMC_set_lp_Yconst_Bernoulli_logit(ERGMM_MCMC_Model *model){
  model->lp_Yconst=0;
}

double ERGMM_MCMC_E_edge_Bernoulli_logit(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j){
  double eta=ERGMM_MCMC_etaij(model,par,i,j);
  return(1/(1+exp(-eta)));
}

/* 
   1 binomial_logit
*/
double ERGMM_MCMC_lp_edge_binomial_logit(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j){
  double eta=ERGMM_MCMC_etaij(model,par,i,j);
  return(model->iY[i][j]*eta-model->iconst[0]*log1p(exp(eta)));
}

void ERGMM_MCMC_set_lp_Yconst_binomial_logit(ERGMM_MCMC_Model *model){
  unsigned int i,j;

  model->lp_Yconst=0;

  if(model->dir){
    for(i=0;i<model->verts;i++)
      for(j=0;j<model->verts;j++)
	if(IS_OBSERVABLE(model->observed_ties,i,j))
	  model->lp_Yconst+=lchoose(model->iconst[0],model->iY[i][j]);
  }
  else{
    for(i=0;i<model->verts;i++)
      for(j=0;j<=i;j++)
	if(IS_OBSERVABLE(model->observed_ties,i,j))
	  model->lp_Yconst+=lchoose(model->iconst[0],model->iY[i][j]);
  }
}

double ERGMM_MCMC_E_edge_binomial_logit(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j){
  double eta=ERGMM_MCMC_etaij(model,par,i,j);
  return(model->iconst[0]/(1+exp(-eta)));
}

/* 2 Poisson_log */
double ERGMM_MCMC_lp_edge_Poisson_log(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j){
  double eta=ERGMM_MCMC_etaij(model,par,i,j);
  return(model->iY[i][j]*eta-exp(eta));
}

void ERGMM_MCMC_set_lp_Yconst_Poisson_log(ERGMM_MCMC_Model *model){
  unsigned int i,j;

  model->lp_Yconst=0;

  if(model->dir){
    for(i=0;i<model->verts;i++)
      for(j=0;j<model->verts;j++)
	if(IS_OBSERVABLE(model->observed_ties,i,j))
	  model->lp_Yconst+=-lgammafn(model->iY[i][j]+1);
  }
  else{
    for(i=0;i<model->verts;i++)
      for(j=0;j<=i;j++)
	if(IS_OBSERVABLE(model->observed_ties,i,j))
	  model->lp_Yconst+=-lgammafn(model->iY[i][j]+1);
  }
}

double ERGMM_MCMC_E_edge_Poisson_log(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j){
  double eta=ERGMM_MCMC_etaij(model,par,i,j);
  return(exp(eta));
}

/*
  3 Bernoulli_cont_logit
*/
double ERGMM_MCMC_lp_edge_Bernoulli_cont_logit(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j){
  double eta=ERGMM_MCMC_etaij(model,par,i,j);
  return(model->dY[i][j]*eta-log1p(exp(eta)));
}

void ERGMM_MCMC_set_lp_Yconst_Bernoulli_cont_logit(ERGMM_MCMC_Model *model){
  model->lp_Yconst=0;
}

double ERGMM_MCMC_E_edge_Bernoulli_cont_logit(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j){
  double eta=ERGMM_MCMC_etaij(model,par,i,j);
  return(1/(1+exp(-eta)));
}

/* 
   4 binomial_cont_logit
*/
double ERGMM_MCMC_lp_edge_binomial_cont_logit(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j){
  double eta=ERGMM_MCMC_etaij(model,par,i,j);
  return(model->dY[i][j]*eta-model->iconst[0]*log1p(exp(eta)));
}

void ERGMM_MCMC_set_lp_Yconst_binomial_cont_logit(ERGMM_MCMC_Model *model){
  unsigned int i,j;

  model->lp_Yconst=0;

  if(model->dir){
    for(i=0;i<model->verts;i++)
      for(j=0;j<model->verts;j++)
	if(IS_OBSERVABLE(model->observed_ties,i,j))
	  model->lp_Yconst+=lchoose(model->iconst[0],model->dY[i][j]);
  }
  else{
    for(i=0;i<model->verts;i++)
      for(j=0;j<=i;j++)
	if(IS_OBSERVABLE(model->observed_ties,i,j))
	  model->lp_Yconst+=lchoose(model->iconst[0],model->dY[i][j]);
  }
}

double ERGMM_MCMC_E_edge_binomial_cont_logit(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j){
  double eta=ERGMM_MCMC_etaij(model,par,i,j);
  return(model->iconst[0]/(1+exp(-eta)));
}

/* 5 Poisson_cont_log */
double ERGMM_MCMC_lp_edge_Poisson_cont_log(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j){
  double eta=ERGMM_MCMC_etaij(model,par,i,j);
  return(model->dY[i][j]*eta-exp(eta));
}

void ERGMM_MCMC_set_lp_Yconst_Poisson_cont_log(ERGMM_MCMC_Model *model){
  unsigned int i,j;

  model->lp_Yconst=0;

  if(model->dir){
    for(i=0;i<model->verts;i++)
      for(j=0;j<model->verts;j++)
	if(IS_OBSERVABLE(model->observed_ties,i,j))
	  model->lp_Yconst+=-lgammafn(model->dY[i][j]+1);
  }
  else{
    for(i=0;i<model->verts;i++)
      for(j=0;j<=i;j++)
	if(IS_OBSERVABLE(model->observed_ties,i,j))
	  model->lp_Yconst+=-lgammafn(model->dY[i][j]+1);
  }
}

double ERGMM_MCMC_E_edge_Poisson_cont_log(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j){
  double eta=ERGMM_MCMC_etaij(model,par,i,j);
  return(exp(eta));
}

/* 6 normal_identity */
double ERGMM_MCMC_lp_edge_normal_identity(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j){
  double diff=model->dY[i][j]-ERGMM_MCMC_etaij(model,par,i,j);
  return(-diff*diff/par->dispersion/2-log(par->dispersion)/2);
}

void ERGMM_MCMC_set_lp_Yconst_normal_identity(ERGMM_MCMC_Model *model){
  unsigned int i,j;

  model->lp_Yconst=0;

  if(model->dir){
    for(i=0;i<model->verts;i++)
      for(j=0;j<model->verts;j++)
	if(IS_OBSERVABLE(model->observed_ties,i,j))
	  model->lp_Yconst+=-M_LN_SQRT_2PI;
  }
  else{
    for(i=0;i<model->verts;i++)
      for(j=0;j<=i;j++)
	if(IS_OBSERVABLE(model->observed_ties,i,j))
	  model->lp_Yconst+=-M_LN_SQRT_2PI;
  }
}

double ERGMM_MCMC_E_edge_normal_identity(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j){
  double eta=ERGMM_MCMC_etaij(model,par,i,j);
  return(eta);
}
