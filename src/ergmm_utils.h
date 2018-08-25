/*  File src/ergmm_utils.h in package latentnet, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2018 Statnet Commons
 */
#ifndef ERGMM_UTILS_H
#define ERGMM_UTILS_H

#include "ergmm_structs.h"
#define rdunif(a,b) ((int) floor(runif(a,b+1)))

double *latentpos_average(double **A, unsigned int n, unsigned int m, double *avA);
void latentpos_translate(double **A, unsigned int n, unsigned int m, double *by);
void randeff_translate(double *v, unsigned int n, double by);
void add_randeff(double *effect, unsigned int n, double **eta, unsigned int is_col);
unsigned int *runifperm(unsigned int n, unsigned int *a);
/*R_INLINE*/ void iswap(int *a, int *b);
/*R_INLINE*/ void uiswap(unsigned int *a, unsigned int *b);
void copy_MCMC_Par(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *source, ERGMM_MCMC_Par *dest);
#endif
