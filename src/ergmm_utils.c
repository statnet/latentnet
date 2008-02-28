/***********************************************************************/
/* Utility functions, mostly involving higher-level matrix operations. */
/***********************************************************************/
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include "matrix_utils.h"
#include "ergmm_utils.h"

// Computes distance between every pair of points and puts the result into dist.
void pairwise_dist(double **A,unsigned int n,unsigned int dim, double **dist){
  unsigned int i,j;
  for(i=0;i<n;i++)
    for(j=0;j<i;j++)
      dist[i][j]=dist[j][i]=dvector_dist(A[i],A[j],dim);
}

// Updates distance between vertex i and all others.
void update_dist(double **A,unsigned int i, unsigned int n, unsigned int dim, double **dist){
  unsigned int j,k;
  double temp,temp2;
  for(j=0;j<n;j++){
    temp2=0;
    for(k=0;k<dim;k++){
      temp=A[i][k]-A[j][k];
      temp2+=temp*temp;
    }
    dist[i][j]=dist[j][i]=sqrt(temp2);
  }
}

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
  if(tocopy(Z)) copy_dmatrix(source->Z,dest->Z,model->verts,model->latent);
  if(tocopy(coef)) copy_dvector(source->coef,dest->coef,model->coef);
  if(tocopy(Z_mean)) copy_dmatrix(source->Z_mean,dest->Z_mean,model->clusters,model->latent);
  if(tocopy(Z_var)) copy_dvector(source->Z_var,dest->Z_var,model->clusters?model->clusters:1);
  if(tocopy(Z_pK)) copy_dvector(source->Z_pK,dest->Z_pK,model->clusters);
  if(tocopy(Z_K)) copy_ivector((int *) source->Z_K,(int *) dest->Z_K,model->verts);
#undef tocopy

  dest->llk=source->llk;
  // The lpedge matrix is NOT copied.
  dest->lpZ=source->lpZ;
  dest->lpLV=source->lpLV;
  dest->lpcoef=source->lpcoef;
}

