#ifndef ERGMM_UTILS_H
#define ERGMM_UTILS_H

#include "ergmm_structs.h"
#define rdunif(a,b) ((int) floor(runif(a,b+1)))

void pairwise_dist(double **A,unsigned int n,unsigned int dim, double **dist);
void update_dist(double **A,unsigned int i, unsigned int n, unsigned int dim, double **dist);
double *latentpos_average(double **A, unsigned int n, unsigned int m, double *avA);
void latentpos_translate(double **A, unsigned int n, unsigned int m, double *by);
unsigned int *runifperm(unsigned int n, unsigned int *a);
/*R_INLINE*/ void iswap(int *a, int *b);
/*R_INLINE*/ void uiswap(unsigned int *a, unsigned int *b);
void copy_MCMC_Par(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *source, ERGMM_MCMC_Par *dest);
#endif
