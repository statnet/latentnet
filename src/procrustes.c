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

  double **Z=dmatrix(*n,*d), **Z_mean=vZ_mean_mcmc?dmatrix(*G,*d):NULL, **Zo=Runpack_dmatrix(vZo,*n,*d,NULL);
  double **A, **tZ, **tZo, **Ahalf, **AhalfInv;
  double **tptrans, **eAvectors, **eADvalues, **teAvectors, *avZ; 
  double *eAvalues, **dd_helper, **dn_helper, **dd2_helper;
  double *workspace;
  
  procr_alloc(*n, *d, *G,
	      &A, &tZ, &tZo, &Ahalf, 
	      &AhalfInv, &tptrans, &eAvectors, 
	      &eADvalues, &teAvectors, &avZ, 
	      &eAvalues, &dd_helper, 
	      &dn_helper, &dd2_helper,
	      &workspace);

  if(*verbose>1) Rprintf("Procrustes: Rotating.\n");
  for(unsigned int s=0; s<*S; s++){
    Runpack_dmatrixs(vZ_mcmc,*n,*d,Z,*S);
    if(vZ_mean_mcmc) Runpack_dmatrixs(vZ_mean_mcmc,*G,*d,Z_mean,*S);
    procr_transform(Z,Z_mean,Zo,*n,*d,*G,
		    Z,Z_mean,
		    A,tZ,tZo,Ahalf,
		    AhalfInv,tptrans,eAvectors,
		    eADvalues,teAvectors,avZ,
		    eAvalues,dd_helper,
		    dn_helper,dd2_helper,workspace);
    Rpack_dmatrixs(Z,*n,*d,vZ_mcmc++,*S);
    if(vZ_mean_mcmc) Rpack_dmatrixs(Z_mean,*G,*d,vZ_mean_mcmc++,*S);
    R_CheckUserInterrupt();
    if(*verbose>2 && (s+1)%(*S/(*verbose))==0) Rprintf("Procrustes: Completed %u/%d.\n",s+1,*S);
  }

  if(*verbose>1) Rprintf("Procrustes: Finished.\n");

  P_free_all();
    
  // Ralloc and friends take care of freeing the memory.
}


void procr_alloc(int n, int d, int G,
		 double ***A, double ***tZ, double ***tZo, double ***Ahalf, 
		 double ***AhalfInv, double ***tptrans, double ***eAvectors, 
		 double ***eADvalues, double ***teAvectors, double **avZ, 
		 double **eAvalues, double ***dd_helper, 
		 double ***dn_helper, double ***dd2_helper,
		 double **workspace){

  // Better allocate too much than too little.
  if(d<G) d=G;

  *dn_helper=dmatrix(d,n);
  *dd_helper=dmatrix(d,d);
  *dd2_helper=dmatrix(d,d);
  *avZ = dvector(d);
  *A=dmatrix(d,d);
  *tZ=dmatrix(d,n);
  *tZo=dmatrix(d,n);
  *eAvectors=dmatrix(d,d);
  *eAvalues=dvector(d);
  *eADvalues=dmatrix(d,d);
  *teAvectors=dmatrix(d,d);
  *Ahalf=*AhalfInv=dmatrix(d,d);
  *tptrans=dmatrix(d,n);
  *workspace=dvector(d*d*3+d*4);
}

/* NOTE: Origins of some of the Procrustes code are unclear. */

int procr_transform(double **Z, double **Z_mean, double **Zo, int n, int d, int G,
		    double **pZ, double **pZ_mean,
		    double **A, double **tZ, double **tZo, double **Ahalf, 
		    double **AhalfInv, double **tptrans, double **eAvectors, 
		    double **eADvalues, double **teAvectors, double *avZ, 
		    double *eAvalues, double **dd_helper, 
		    double **dn_helper, double **dd2_helper,
		    double *workspace)
{
  
  int i=0,j=0;
  /*   double *tempZ, *tempZo, *temp ; */

  /*  First centers Z around the origin  */
  
  for(j=0;j<d;j++){
    avZ[j] = 0.0;
    for(i=0;i<n;i++){
      avZ[j] += Z[i][j] / n;
    }
  }

  /*   subtract the averages */ 
  for(j=0;j<d;j++){
    for(i=0;i<n;i++){
      Z[i][j] -= avZ[j];
    }
  }

  /* Compute A = tZ*Zo*tZo*Z*/
  t(Z,n,d,tZ);
  t(Zo,n,d,tZo);

  init_dmatrix(A,d,d,0.0);
  init_dmatrix(dd_helper,d,d,0.0);
  init_dmatrix(dn_helper,d,n,0.0);
  dmatrix_multiply(tZ,d,n,Zo,d,dd_helper);
  dmatrix_multiply(dd_helper,d,d,tZo,n,dn_helper);
  dmatrix_multiply(dn_helper,d,n,Z,d,A);

  /* Compute sqrt(A) */
  sym_eigen(A,d,1,eAvalues,eAvectors,workspace);
  for(i=0;i<d;i++){
    eADvalues[i][i] = sqrt(eAvalues[i]);
  }
  t(eAvectors,d,d,teAvectors);  

  init_dmatrix(dd_helper,d,d,0.0);
  init_dmatrix(Ahalf,d,d,0.0);
  dmatrix_multiply(eAvectors,d,d,eADvalues,d,dd_helper);
  dmatrix_multiply(dd_helper,d,d,teAvectors,d,Ahalf);

  /*   init_dmatrix(AhalfInv,d,d,0.0); */
  /* Now compute the inverse */
  if(inverse(Ahalf,d,AhalfInv,workspace) != FOUND) return NOTFOUND;

  /*Now compute t(Zo)*Z*AhalfInv*tZ*/
  init_dmatrix(dd_helper,d,d,0.0);
  init_dmatrix(dd2_helper,d,d,0.0);
  init_dmatrix(tptrans,d,n,0.0);
  dmatrix_multiply(tZo,d,n,Z,d,dd_helper);
  dmatrix_multiply(dd_helper,d,d,AhalfInv,d,dd2_helper);

  // Add the center back on to tZ.
  for(j=0;j<d;j++){
    for(i=0;i<n;i++){
      tZ[j][i] += avZ[j];
    }
  }

  dmatrix_multiply(dd2_helper,d,d,tZ,n,tptrans);
  t(tptrans,d,n,pZ);

  if(Z_mean){
    // Note that n>=G, so we can reuse tZ et al for Z_mean.
    t(Z_mean,G,d,tZ);
    init_dmatrix(tptrans,d,n,0.0);
    dmatrix_multiply(dd2_helper,d,d,tZ,G,tptrans);
    t(tptrans,d,G,pZ_mean);
  }

  return(FOUND);
}
