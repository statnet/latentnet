/*  File src/ergmm_latent_effects.h in package latentnet, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2018 Statnet Commons
 */
#ifndef ERGMM_LATENT_EFFECTS_H
#define ERGMM_LATENT_EFFECTS_H

#define N_LATENT_EFF 3

/* Declare "lookup tables" for latent effects. */
typedef double (*ERGMM_MCMC_latent_eff_t)(double *u, double *v, unsigned int dim);

extern ERGMM_MCMC_latent_eff_t ERGMM_MCMC_latent_eff[];

/* Latent effect # */
/* 0 */ 
double ERGMM_MCMC_latent_eff_negative_Euclidean_distance(double *u, double *v, unsigned int dim);
/* 1 */ 
double ERGMM_MCMC_latent_eff_dot_product(double *u, double *v, unsigned int dim);
/* 2 */ 
double ERGMM_MCMC_latent_eff_negative_Euclidean_distance2(double *u, double *v, unsigned int dim);
#endif/* ERGMM_LATENT_EFFECTS_H */
