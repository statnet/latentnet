/*  File src/procrustes.c in package latentnet, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2014 Statnet Commons
 */
/********************************************/
/* Utilities for processing posterior draws */
/********************************************/

#include <R.h>
#include <Rmath.h>
#include "procrustes.h"
#include "matrix_utils.h"

void procr_transform_wrapper(int *S, int *n, int *d, int *G, double *vZo,
			     double *vZ_mcmc, double *vZ_mean_mcmc, int *verbose){

  if(*verbose>1) Rprintf("Procrustes: Allocating memory.\n");

  double **Z=dmatrix(*n,*d), **pZ=dmatrix(*n,*d), **Z_mean=*G>0?dmatrix(*G,*d):NULL, **pZ_mean=*G>0?dmatrix(*G,*d):NULL, **Zo=Runpack_dmatrix(vZo,*n,*d,NULL);
  double *avZ, **tZZo, **U, *Sigma, **tV, *workspace;
  
  procr_alloc(*n, *d, *G,
	      &avZ, &tZZo, &U, &Sigma, &tV, &workspace);

  if(*verbose>1) Rprintf("Procrustes: Rotating.\n");
  for(unsigned int s=0; s<*S; s++){
    Runpack_dmatrixs(vZ_mcmc,*n,*d,Z,*S);
    if(vZ_mean_mcmc) Runpack_dmatrixs(vZ_mean_mcmc,*G,*d,Z_mean,*S);
    procr_transform(Z,Z_mean,Zo,*n,*d,*G,
		    pZ,pZ_mean, avZ, tZZo, U, Sigma, tV, workspace);
    Rpack_dmatrixs(pZ,*n,*d,vZ_mcmc++,*S);
    if(vZ_mean_mcmc) Rpack_dmatrixs(pZ_mean,*G,*d,vZ_mean_mcmc++,*S);
    R_CheckUserInterrupt();
    if(*verbose>2 && (s+1)%(*S/(*verbose))==0) Rprintf("Procrustes: Completed %u/%d.\n",s+1,*S);
  }

  if(*verbose>1) Rprintf("Procrustes: Finished.\n");

  P_free_all();
    
  // Ralloc and friends take care of freeing the memory.
}


void procr_alloc(int n, int d, int G,
		 double **avZ, double ***tZZo, double ***U, double **S, double ***tV,
		 double **workspace){

  // Better allocate too much than too little.
  if(d<G) d=G;

  *avZ = dvector(d);
  *tZZo=dmatrix(d,d);
  *U=dmatrix(d,d);
  *S=dvector(d);
  *tV=dmatrix(d,d);
  *workspace=dvector(d*d*5+d*5);
}

/* NOTE: Origins of some of the Procrustes code are unclear. */

int procr_transform(double **Z, double **Z_mean, double **Zo, int n, int d, int G,
		    double **pZ, double **pZ_mean,
		    double *avZ, double **tZZo, double **U, double *S, double **tV,
		    double *workspace)
{
  
  /*  First center Z around the origin  */
  
  for(unsigned int j=0;j<d;j++){
    avZ[j] = 0.0;
    for(unsigned int i=0;i<n;i++){
      avZ[j] += Z[i][j] / n;
    }
  }

  /*   subtract the averages */ 
  for(unsigned int j=0;j<d;j++){
    for(unsigned int i=0;i<n;i++){
      Z[i][j] -= avZ[j];
    }
  }

  /* Compute  tZ*Zo */
  dmatrix_init(tZZo, d, d, 0.0);
  dmatrix_crossprod(Z,n,d,Zo,d,tZZo);

  if(dgesvd_full_wrapper(tZZo, d, d, U, S, tV, workspace, d*d*5+d*5)!=0) return NOTFOUND;

  double **R = tZZo; // Reuse storage for tZZo for the rotation matrix.

  dmatrix_init(R, d, d, 0.0);
  dmatrix_multiply(U,d,d,tV,d,R);

  // Add the center back on to tZ.
  for(unsigned int j=0;j<d;j++){
    for(unsigned int i=0;i<n;i++){
      Z[i][j] += avZ[j];
    }
  }

  dmatrix_init(pZ, n, d, 0.0);
  dmatrix_multiply(Z,n,d,R,d,pZ);

  if(Z_mean){
    dmatrix_init(pZ_mean, G, d, 0.0);
    dmatrix_multiply(Z_mean,G,d,R,d, pZ_mean);
  }
  
  return(FOUND);
}
