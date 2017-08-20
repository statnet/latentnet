/*  File src/ergmm_latent_effects.c in package latentnet, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2017 Statnet Commons
 */
/********************************************************/
/* Types of latent space effects supported by latentnet */
/********************************************************/
#include "ergmm_latent_effects.h"
#include <math.h>
/* ERGMM_MCMC_latent_eff = effect of the latent space terms on the linear component. */

/* Define "lookup tables" for latent space effects. */

double (*ERGMM_MCMC_latent_eff[N_LATENT_EFF])(double *u, double *v, unsigned int dim)={
  ERGMM_MCMC_latent_eff_negative_Euclidean_distance,
  ERGMM_MCMC_latent_eff_dot_product,
  ERGMM_MCMC_latent_eff_negative_Euclidean_distance2
};

/* 
   0 Negative Euclidean distance
*/ 
double ERGMM_MCMC_latent_eff_negative_Euclidean_distance(double *u, double *v, unsigned int dim){
  unsigned int k;
  double dist,dist2=0;
  for(k=0;k<dim;k++){
    dist=u[k]-v[k];
    dist2+=dist*dist;
  }
  return(-sqrt(dist2));
}


/*
  1 Dot product (bilinear) effect
*/
double ERGMM_MCMC_latent_eff_dot_product(double *u, double *v, unsigned int dim){
  double prod=0;
  for(unsigned int i=0; i<dim; i++)
    prod += u[i]*v[i];
  return(prod);
}

/* 
   2 Negative Euclidean distance squared
*/ 
double ERGMM_MCMC_latent_eff_negative_Euclidean_distance2(double *u, double *v, unsigned int dim){
  unsigned int k;
  double dist,dist2=0;
  for(k=0;k<dim;k++){
    dist=u[k]-v[k];
    dist2+=dist*dist;
  }
  return(-dist2);
}
