#ifndef ERGMM_LATENT_EFFECTS_H
#define ERGMM_LATENT_EFFECTS_H

#define N_LATENT_EFF 3

/* Declare "lookup tables" for latent effects. */
double (*ERGMM_MCMC_latent_eff[N_LATENT_EFF])(double *u, double *v, unsigned int dim);

/* Latent effect # */
/* 0 */ 
double ERGMM_MCMC_latent_eff_negative_Euclidean_distance(double *u, double *v, unsigned int dim);
/* 1 */ 
double ERGMM_MCMC_latent_eff_dot_product(double *u, double *v, unsigned int dim);
/* 2 */ 
double ERGMM_MCMC_latent_eff_negative_Euclidean_distance2(double *u, double *v, unsigned int dim);
#endif/* ERGMM_LATENT_EFFECTS_H */
