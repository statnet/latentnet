#include "ergmm_families.h"
#include "ergmm_probs.h"
#include <R.h>
#include <Rmath.h>
#include <math.h>

/*
  Bernoulli-logit
*/
double ERGMM_MCMC_lp_edge_Bernoulli_logit(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j){
  double eta=ERGMM_MCMC_etaij(model,par,i,j);
  return(par->lpedge[i][j]=model->Y[i][j]*eta-log(1+exp(eta)));
}

void ERGMM_MCMC_set_lp_Yconst_Bernoulli_logit(ERGMM_MCMC_Model *model){
  model->lp_Yconst=0;
}

double ERGMM_MCMC_E_edge_Bernoulli_logit(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j){
  double eta=ERGMM_MCMC_etaij(model,par,i,j);
  return(1/(1+exp(-eta)));
}

/* 
   binomial-logit
*/
double ERGMM_MCMC_lp_edge_binomial_logit(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j){
  double eta=ERGMM_MCMC_etaij(model,par,i,j);
  return(par->lpedge[i][j]=model->Y[i][j]*eta-model->iconst[0]*log(1+exp(eta)));
}

void ERGMM_MCMC_set_lp_Yconst_binomial_logit(ERGMM_MCMC_Model *model){
  unsigned int i,j;

  model->lp_Yconst=0;

  if(model->dir){
    for(i=0;i<model->verts;i++)
      for(j=0;j<model->verts;j++)
	if(IS_OBSERVABLE(model->observed_ties,i,j))
	  model->lp_Yconst+=lchoose(model->iconst[0],model->Y[i][j]);
  }
  else{
    for(i=0;i<model->verts;i++)
      for(j=0;j<=i;j++)
	if(IS_OBSERVABLE(model->observed_ties,i,j))
	  model->lp_Yconst+=lchoose(model->iconst[0],model->Y[i][j]);
  }
}

double ERGMM_MCMC_E_edge_binomial_logit(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j){
  double eta=ERGMM_MCMC_etaij(model,par,i,j);
  return(model->iconst[0]/(1+exp(-eta)));
}

/* Poisson-log */
double ERGMM_MCMC_lp_edge_Poisson_log(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j){
  double eta=ERGMM_MCMC_etaij(model,par,i,j);
  return(par->lpedge[i][j]=model->Y[i][j]*eta-exp(eta));
}

void ERGMM_MCMC_set_lp_Yconst_Poisson_log(ERGMM_MCMC_Model *model){
  unsigned int i,j;

  model->lp_Yconst=0;

  if(model->dir){
    for(i=0;i<model->verts;i++)
      for(j=0;j<model->verts;j++)
	if(IS_OBSERVABLE(model->observed_ties,i,j))
	  model->lp_Yconst+=-lgammafn(model->Y[i][j]+1);
  }
  else{
    for(i=0;i<model->verts;i++)
      for(j=0;j<=i;j++)
	if(IS_OBSERVABLE(model->observed_ties,i,j))
	  model->lp_Yconst+=-lgammafn(model->Y[i][j]+1);
  }
}

double ERGMM_MCMC_E_edge_Poisson_log(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j){
  double eta=ERGMM_MCMC_etaij(model,par,i,j);
  return(exp(eta));
}
