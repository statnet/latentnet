/*  File src/ergmm_updaters.h in package latentnet, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2020 Statnet Commons
 */
#ifndef ERGMM_UPDATERS_H
#define ERGMM_UPDATERS_H

#include "ergmm_structs.h"
void ERGMM_MCMC_propose(ERGMM_MCMC_Model *model, ERGMM_MCMC_MCMCState *cur, unsigned int Z, unsigned int RE, unsigned int coef, unsigned int LV, unsigned int REV, unsigned int dispersion);
/* The main difference between accepting and rejecting is the direction, so
   they are just macros. (The other main difference is whether lpedge needs
   to be copied.)
 */
#define ERGMM_MCMC_accept(model,cur) ERGMM_MCMC_prop_end(model,cur,cur->prop,cur->state, TRUE)
#define ERGMM_MCMC_reject(model,cur) ERGMM_MCMC_prop_end(model,cur,cur->state,cur->prop, FALSE)
void ERGMM_MCMC_prop_end(ERGMM_MCMC_Model *model, ERGMM_MCMC_MCMCState *cur, ERGMM_MCMC_Par *new, ERGMM_MCMC_Par *old, unsigned int copy_lpedge);
unsigned int ERGMM_MCMC_Z_RE_up(ERGMM_MCMC_Model *model, ERGMM_MCMC_Priors *prior, ERGMM_MCMC_MCMCState *cur,ERGMM_MCMC_MCMCSettings *setting);
unsigned int ERGMM_MCMC_coef_up_scl_Z_shift_RE(ERGMM_MCMC_Model *model,  ERGMM_MCMC_Priors *prior, ERGMM_MCMC_MCMCState *cur,ERGMM_MCMC_MCMCSettings *setting);
void ERGMM_MCMC_CV_up(ERGMM_MCMC_Model *model, ERGMM_MCMC_Priors *prior, ERGMM_MCMC_MCMCState *cur);
void ERGMM_MCMC_LV_up(ERGMM_MCMC_Model *model, ERGMM_MCMC_Priors *prior, ERGMM_MCMC_MCMCState *cur);
void ERGMM_MCMC_REV_up(ERGMM_MCMC_Model *model, ERGMM_MCMC_Priors *prior, ERGMM_MCMC_MCMCState *cur);

#endif
