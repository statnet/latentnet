#include "ergmm_families.h"
#include "ergmm_probs.h"
#include <R.h>
#include <Rmath.h>
#include <math.h>


/* Define "lookup tables" for families. */

const unsigned int ERGMM_MCMC_is_discrete[6]={TRUE,TRUE,TRUE,FALSE,FALSE,FALSE};
const unsigned int ERGMM_MCMC_to_cont[6]={3,4,5,3,4,5};

double (*ERGMM_MCMC_lp_edge[N_FAMILIES])(ERGMM_MCMC_Model *, ERGMM_MCMC_Par *,
					 unsigned int, unsigned int)={
  ERGMM_MCMC_lp_edge_Bernoulli_logit,
  ERGMM_MCMC_lp_edge_binomial_logit,
  ERGMM_MCMC_lp_edge_Poisson_log,
  ERGMM_MCMC_lp_edge_Bernoulli_cont_logit,
  ERGMM_MCMC_lp_edge_binomial_cont_logit,
  ERGMM_MCMC_lp_edge_Poisson_cont_log
};
  
void (*ERGMM_MCMC_set_lp_Yconst[N_FAMILIES])(ERGMM_MCMC_Model *)={
  ERGMM_MCMC_set_lp_Yconst_Bernoulli_logit,
  ERGMM_MCMC_set_lp_Yconst_binomial_logit,
  ERGMM_MCMC_set_lp_Yconst_Poisson_log,
  ERGMM_MCMC_set_lp_Yconst_Bernoulli_cont_logit,
  ERGMM_MCMC_set_lp_Yconst_binomial_cont_logit,
  ERGMM_MCMC_set_lp_Yconst_Poisson_cont_log
};
  
double (*ERGMM_MCMC_E_edge[N_FAMILIES])(ERGMM_MCMC_Model *, ERGMM_MCMC_Par *,
					unsigned int, unsigned int)={
  ERGMM_MCMC_E_edge_Bernoulli_logit,
  ERGMM_MCMC_E_edge_binomial_logit,
  ERGMM_MCMC_E_edge_Poisson_log,
  ERGMM_MCMC_E_edge_Bernoulli_cont_logit,
  ERGMM_MCMC_E_edge_binomial_cont_logit,
  ERGMM_MCMC_E_edge_Poisson_cont_log
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
  return(model->dY[i][j]*eta-model->dconst[0]*log1p(exp(eta)));
}

void ERGMM_MCMC_set_lp_Yconst_binomial_cont_logit(ERGMM_MCMC_Model *model){
  unsigned int i,j;

  model->lp_Yconst=0;

  if(model->dir){
    for(i=0;i<model->verts;i++)
      for(j=0;j<model->verts;j++)
	if(IS_OBSERVABLE(model->observed_ties,i,j))
	  model->lp_Yconst+=lchoose(model->dconst[0],model->dY[i][j]);
  }
  else{
    for(i=0;i<model->verts;i++)
      for(j=0;j<=i;j++)
	if(IS_OBSERVABLE(model->observed_ties,i,j))
	  model->lp_Yconst+=lchoose(model->dconst[0],model->dY[i][j]);
  }
}

double ERGMM_MCMC_E_edge_binomial_cont_logit(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j){
  double eta=ERGMM_MCMC_etaij(model,par,i,j);
  return(model->dconst[0]/(1+exp(-eta)));
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
