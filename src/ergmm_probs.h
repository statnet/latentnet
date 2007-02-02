#ifndef LATENT_PROBS_H
#define LATENT_PROBS_H 1
#include <R.h>

#include "ergmm_structs.h"

#define LOG_ILOGIT(Y,eta,i,j) (Y[i][j]*eta[i][j]-log(1+exp(eta[i][j])))
#define IS_OBSERVABLE(obs_ties,i,j) (obs_ties ? obs_ties[i][j] : i!=j)

R_INLINE double ERGMM_MCMC_etaij(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par,unsigned int i,unsigned int j);
double ERGMM_MCMC_lp_Y(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int own_lpedge);
double ERGMM_MCMC_lp_Y_diff(ERGMM_MCMC_Model *model, ERGMM_MCMC_MCMCState *cur);
double ERGMM_MCMC_logp_Z(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par);
double ERGMM_MCMC_logp_Z_diff(ERGMM_MCMC_Model *model, ERGMM_MCMC_MCMCState *cur);
double ERGMM_MCMC_logp_RE(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par);
double ERGMM_MCMC_logp_RE_diff(ERGMM_MCMC_Model *model, ERGMM_MCMC_MCMCState *cur);
double ERGMM_MCMC_logp_LV(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, ERGMM_MCMC_Priors *prior);
double ERGMM_MCMC_logp_LV_diff(ERGMM_MCMC_Model *model, ERGMM_MCMC_MCMCState *cur, ERGMM_MCMC_Priors *prior);
double ERGMM_MCMC_logp_coef(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, ERGMM_MCMC_Priors *prior);
double ERGMM_MCMC_logp_coef_diff(ERGMM_MCMC_Model *model, ERGMM_MCMC_MCMCState *cur, ERGMM_MCMC_Priors *prior);
double ERGMM_MCMC_logp_REV(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, ERGMM_MCMC_Priors *prior);
double ERGMM_MCMC_logp_REV_diff(ERGMM_MCMC_Model *model, ERGMM_MCMC_MCMCState *cur, ERGMM_MCMC_Priors *prior);

void ERGMM_lp_Y_wrapper(int *n, int *p, int *d,
			   int *dir, int *vY,
			   int *family, int *iconsts, double *dconsts,
			   double *vX, double *vZ,
			   double *coef,
			   double *sender, double *receiver, int *lock_RE,
			   int *vobserved_ties,
			   double *llk);


#endif
