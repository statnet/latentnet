/*  File src/ergmm_utils.c in package latentnet, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2022 Statnet Commons
 */
/***********************************************************************/
/* Utility functions, mostly involving higher-level matrix operations. */
/***********************************************************************/
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include "matrix_utils.h"
#include "ergmm_utils.h"

double *latentpos_average(double **A, unsigned int n, unsigned int m, double *avA){
  unsigned int i,j;
  if(!avA) avA=dvector(m);
  init_dvector(avA,m,0);
  for(j=0;j<m;j++){
    for(i=0;i<n;i++)
      avA[j]+=A[i][j];
    avA[j]/=n;
  }
  return(avA);
}

void latentpos_translate(double **A, unsigned int n, unsigned int m, double *by){
  unsigned int i,j;
  for(j=0;j<m;j++)
    for(i=0;i<n;i++)
      A[i][j]+=by[j];
}

void randeff_translate(double *v, unsigned int n, double by){
  unsigned int i;
  for(i=0;i<n;i++) v[i]+=by;
}

void add_randeff(double *effect, unsigned int n, double **eta, unsigned int is_col){
  unsigned int i,j;
  if(is_col)
    for(i=0;i<n;i++){
      for(j=0;j<n;j++){
	eta[i][j]+=effect[i];
      }
    }
  else
    for(i=0;i<n;i++){
      for(j=0;j<n;j++){
	eta[i][j]+=effect[j];
      }
    }
}

/* Generate a uniformly random permutation. */
unsigned int *runifperm(unsigned int n, unsigned int *a){
  unsigned int i;
  if(!a) a=(unsigned int *) ivector(n);
  
  for(i=0;i<n;i++) a[i]=i;

  for(i=0;i<n-1;i++) uiswap(a+i, a+rdunif(i,n-1));

  return(a);
}

/*R_INLINE*/ void iswap(int *a, int *b){
  int tmp=*b;
  *b=*a;
  *a=tmp;
}

/*R_INLINE*/ void uiswap(unsigned int *a, unsigned int *b){
  unsigned int tmp=*b;
  *b=*a;
  *a=tmp;
}

void copy_MCMC_Par(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *source, ERGMM_MCMC_Par *dest){
#define tocopy(name) (source->name && (source->name != dest->name))
  if(tocopy(Z)) dmatrix_copy_contents(source->Z,dest->Z,model->verts,model->latent);
  if(tocopy(coef)) copy_dvector(source->coef,dest->coef,model->coef);
  if(tocopy(Z_mean)) dmatrix_copy_contents(source->Z_mean,dest->Z_mean,model->clusters,model->latent);
  if(tocopy(Z_var)) copy_dvector(source->Z_var,dest->Z_var,model->clusters?model->clusters:1);
  if(tocopy(Z_pK)) copy_dvector(source->Z_pK,dest->Z_pK,model->clusters);
  if(tocopy(sender)) copy_dvector(source->sender,dest->sender,model->verts);
  if(source->sender) dest->sender_var=source->sender_var;
  if(!model->sociality && tocopy(receiver)) copy_dvector(source->receiver,dest->receiver,model->verts);
  if(source->receiver) dest->receiver_var=source->receiver_var;
  if(model->dispersion) dest->dispersion=source->dispersion;
  if(tocopy(Z_K)) copy_ivector((int *) source->Z_K,(int *) dest->Z_K,model->verts);
#undef tocopy

  dest->llk=source->llk;
  // The lpedge matrix is NOT copied.
  dest->lpZ=source->lpZ;
  dest->lpLV=source->lpLV;
  dest->lpcoef=source->lpcoef;
  dest->lpRE=source->lpRE;
  dest->lpREV=source->lpREV;
}

