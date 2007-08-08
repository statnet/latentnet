#ifndef ERGMM_FAMILIES_H
#define ERGMM_FAMILIES_H

#include "ergmm_structs.h"

#ifndef TRUE
#define TRUE !0
#endif /* TRUE */

#ifndef FALSE
#define FALSE 0
#endif /* FALSE */

#define N_FAMILIES 6
/* Declare "lookup tables" for families. */

const unsigned int ERGMM_MCMC_is_discrete[N_FAMILIES];
const unsigned int ERGMM_MCMC_to_cont[N_FAMILIES];

double (*ERGMM_MCMC_lp_edge[N_FAMILIES])(ERGMM_MCMC_Model *, ERGMM_MCMC_Par *,
					 unsigned int, unsigned int);
void (*ERGMM_MCMC_set_lp_Yconst[N_FAMILIES])(ERGMM_MCMC_Model *);  
double (*ERGMM_MCMC_E_edge[N_FAMILIES])(ERGMM_MCMC_Model *, ERGMM_MCMC_Par *,
					unsigned int, unsigned int);

/* Family # */
/* 0 */
double ERGMM_MCMC_lp_edge_Bernoulli_logit(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j);
void ERGMM_MCMC_set_lp_Yconst_Bernoulli_logit(ERGMM_MCMC_Model *model);
double ERGMM_MCMC_E_edge_Bernoulli_logit(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j);

/* 1 */
double ERGMM_MCMC_lp_edge_binomial_logit(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j);
void ERGMM_MCMC_set_lp_Yconst_binomial_logit(ERGMM_MCMC_Model *model);
double ERGMM_MCMC_E_edge_binomial_logit(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j);

/* 2 */
double ERGMM_MCMC_lp_edge_Poisson_log(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j);
void ERGMM_MCMC_set_lp_Yconst_Poisson_log(ERGMM_MCMC_Model *model);
double ERGMM_MCMC_E_edge_Poisson_log(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j);

/* 3 */
double ERGMM_MCMC_lp_edge_Bernoulli_cont_logit(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j);
void ERGMM_MCMC_set_lp_Yconst_Bernoulli_cont_logit(ERGMM_MCMC_Model *model);
double ERGMM_MCMC_E_edge_Bernoulli_cont_logit(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j);

/* 4 */
double ERGMM_MCMC_lp_edge_binomial_cont_logit(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j);
void ERGMM_MCMC_set_lp_Yconst_binomial_cont_logit(ERGMM_MCMC_Model *model);
double ERGMM_MCMC_E_edge_binomial_cont_logit(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j);

/* 5 */
double ERGMM_MCMC_lp_edge_Poisson_cont_log(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j);
void ERGMM_MCMC_set_lp_Yconst_Poisson_cont_log(ERGMM_MCMC_Model *model);
double ERGMM_MCMC_E_edge_Poisson_cont_log(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par, unsigned int i, unsigned int j);

#endif /* ERGMM_FAMILIES_H */
