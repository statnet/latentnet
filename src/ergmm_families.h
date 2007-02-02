#ifndef ERGMM_FAMILIES_H
#define ERGMM_FAMILIES_H 1

#include "ergmm_structs.h"

/* Family # */
/* 0 */
double ERGMM_MCMC_lp_edge_Bernoulli_logit(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j);
void ERGMM_MCMC_set_lp_Yconst_Bernoulli_logit(ERGMM_MCMC_Model *model);
/* 1 */
double ERGMM_MCMC_lp_edge_binomial_logit(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j);
void ERGMM_MCMC_set_lp_Yconst_binomial_logit(ERGMM_MCMC_Model *model);
/* 2 */
double ERGMM_MCMC_lp_edge_Poisson_log(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j);
void ERGMM_MCMC_set_lp_Yconst_Poisson_log(ERGMM_MCMC_Model *model);


#endif
