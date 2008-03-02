/****************/
/* KL Switching */
/****************/

#include <R.h>
#include <Rmath.h>
#include "kl_switch.h"
#include "matrix_utils.h"
#include "mvnorm.h"

void klswitch_wrapper(int *maxit, int *S, int *n, int *d, int *G,
		      double *vZ_mcmc, int *Z_ref, double *vZ_mean_mcmc, double *vZ_var_mcmc,
		      int *vZ_K_mcmc, double *vZ_pK_mcmc,
		      double *vQ, int *verbose){

  /* R function enabling uniform RNG (for tie-beaking) */
  GetRNGstate();
  

  if(*verbose>1) Rprintf("KLswitch: Allocating memory.\n");
  ERGMM_MCMC_Par *sample = (ERGMM_MCMC_Par *) P_alloc(*S,sizeof(ERGMM_MCMC_Par));
  unsigned int *perm=(unsigned int *)P_alloc(*G,sizeof(unsigned int)), *bestperm=(unsigned int *)P_alloc(*G,sizeof(unsigned int));
  unsigned int *dir=(unsigned int *)P_alloc(*G,sizeof(unsigned int));
  double **Q=Runpack_dmatrix(vQ,*n,*G,NULL);
  ERGMM_MCMC_Par tmp;
    
  tmp.Z_mean=dmatrix(*G,*d);
  tmp.Z_var=dvector(*G);
  tmp.Z_pK=dvector(*G);
  tmp.Z_K=ivector(*n);

  double  ***pK=d3array(*S,*n,*G), 
    ***Z_space=*Z_ref ? d3array(1,*n,*d):d3array(*S,*n,*d),
    ***Z_mean_space=d3array(*S,*G,*d), **Z_var_space=dmatrix(*S,*G),
    **Z_pK_space = dmatrix(*S, *G);
  unsigned int **Z_K_space= (unsigned int **) imatrix(*S,*n);

  if(*Z_ref) Runpack_dmatrix(vZ_mcmc, *n, *d, Z_space[0]);

  if(*verbose>1) Rprintf("KLswitch: Unpacking R input and precalculating pK.\n");

  for(unsigned int s=0; s<*S; s++){
    ERGMM_MCMC_Par *cur=sample+s;
    if(*Z_ref) cur->Z = Z_space[0];
    else cur->Z = Runpack_dmatrixs(vZ_mcmc+s,*n,*d,Z_space[s],*S);
    cur->Z_mean = Runpack_dmatrixs(vZ_mean_mcmc+s,*G,*d,Z_mean_space[s],*S);
    cur->Z_var = Runpack_dvectors(vZ_var_mcmc+s,*G,Z_var_space[s],*S);
    cur->Z_pK = Runpack_dvectors(vZ_pK_mcmc+s,*G,Z_pK_space[s],*S);
    cur->Z_K = Runpack_ivectors(vZ_K_mcmc+s,*n,Z_K_space[s],*S);

    /* Precalculate p(K|Z,mu,sigma).
       This is a speed-memory trade-off:
       we allocate approx. S*n*G*8 bytes so that we don't have to do
       O(it*S*n*G*G!) (!=factorial) evaluations of MVN density.
     */

    for(unsigned int i=0; i<*n; i++){
      double pKsum=0;
      for(unsigned int g=0; g<*G; g++){
	pK[s][i][g]=dindnormmu(*d,cur->Z[i],
			       cur->Z_mean[g],sqrt(cur->Z_var[g]),FALSE);
	pK[s][i][g]*=cur->Z_pK[g];
	pKsum+=pK[s][i][g];
      }
      for(unsigned int g=0; g<*G; g++){
	pK[s][i][g]/=pKsum;
      }
    }
    if(*verbose>2 && (s+1)%(*S/(*verbose))==0) Rprintf("KLswitch: Unpacking & precalculating: Completed %u/%d.\n",s+1,*S);
  }

  if(*verbose>1) Rprintf("KLswitch: Iterating between label-switching to Q and recalculating Q.\n");

  for(unsigned int it=0; it<*maxit; it++){
    /* Do NOT switch the function call and "it>0", or you will short-circuit the function call.

       The first time through the loop, label-switch to MKL of Q, and recalculate Q.
       During subsequent iterations, if no labels were permuted, Q need not be recalculated (
       since it already has been calculated with those labels). 
    */
    if(!klswitch_step2(Q, sample, &tmp, *S, *n, *d, *G, perm, bestperm, dir, pK) && it>0){
      if(*verbose>1) Rprintf("KLswitch: Iterating: Q converged after %u iterations.\n", it+1);
      break;
    }
    klswitch_step1(sample, *S, *n, *d, *G, Q, pK);


    if(*verbose>2) Rprintf("KLswitch: Iterating: Completed %u/%d.\n",it+1,*maxit);
  }

  if(*verbose>1) Rprintf("KLswitch: Packing for return to R.\n");

  for(unsigned int s=0; s<*S; s++){
    ERGMM_MCMC_Par *cur=sample+s;
    Rpack_dmatrixs(cur->Z_mean,*G,*d,vZ_mean_mcmc+s,*S);
    Rpack_dvectors(cur->Z_var,*G,vZ_var_mcmc+s,*S);
    Rpack_dvectors(cur->Z_pK,*G,vZ_pK_mcmc+s,*S);
    Rpack_ivectors(cur->Z_K,*n,vZ_K_mcmc+s,*S);
  }

  Rpack_dmatrixs(Q,*n,*G,vQ,1);

  if(*verbose>1) Rprintf("KLswitch: Finished.\n");

  PutRNGstate();
  P_free_all();

}

// Generate the next permutation using Johnson-Trotter Algorithm.
/*R_INLINE*/ int nextperm(unsigned int n, unsigned int *p, unsigned int *dir){
  unsigned int tmpd;
  unsigned int tmpp=0;
  unsigned int imax=0;

  // dir==0 -> left, dir!=0 -> right.

  // Find the greater "mobile" integer... 
  for(unsigned int i=0;i<n;i++){
    if(((i<n-1 && dir[i] && p[i]>p[i+1])  //       value in i'th position is "mobile" to the right
	||                                //    or
	(i>0 && !dir[i] && p[i]>p[i-1]))  //       value in i'th position is "mobile" to the left
       &&                                 // and
       p[i]>tmpp                          //    it's greater than the current greater "mobile"
       ){
      tmpp=p[i];
      imax=i;
    }
  }
  // If no "mobile" integers, we are done.
  if(!tmpp) return 0;

  tmpd=dir[imax];
  tmpp=p[imax];
  if(tmpd){
    // Swap to the right...
    dir[imax]=dir[imax+1];
    p[imax]=p[imax+1];
    dir[imax+1]=tmpd;
    p[imax+1]=tmpp;
  }else{
    // Swap to the left...
    dir[imax]=dir[imax-1];
    p[imax]=p[imax-1];
    dir[imax-1]=tmpd;
    p[imax-1]=tmpp;
  }
  // Switch direction of those greater than the greatest "mobile".
  for(unsigned int i=0;i<n;i++)
    if(p[i]>tmpp) dir[i]=!dir[i];

  return !0;
} 


/*R_INLINE*/ void apply_perm(unsigned int *perm, ERGMM_MCMC_Par *to, double **pK, ERGMM_MCMC_Par *tmp, int n, int d, int G){
  copy_dmatrix(to->Z_mean,tmp->Z_mean,G,d);
  copy_dvector(to->Z_var,tmp->Z_var,G);
  copy_dvector(to->Z_pK,tmp->Z_pK,G);
  copy_ivector(to->Z_K,tmp->Z_K,n);

  for(unsigned int g=0; g<G; g++){
    copy_dvector(tmp->Z_mean[perm[g]-1],to->Z_mean[g],d);
    to->Z_var[g]=tmp->Z_var[perm[g]-1];
    to->Z_pK[g]=tmp->Z_pK[perm[g]-1];
    for(unsigned int i=0; i<n; i++){
      if(tmp->Z_K[i]==perm[g]) to->Z_K[i]=g+1;
    }
  }

  for(unsigned int i=0; i<n; i++){
    // Reuse tmp->Z_pK to store the assignment probabilities for permutation.
    double *pKtmp=tmp->Z_pK;
    copy_dvector(pK[i],pKtmp,G);
    for(unsigned int g=0; g<G; g++){
      pK[i][g]=pKtmp[perm[g]-1];
    }
  }
}

// "Step 1": evaluate new Q.
void klswitch_step1(ERGMM_MCMC_Par *sample, int S, int n, int d, int G, double **Q, double ***pK){
  for(unsigned int i=0; i<n; i++){
    for(unsigned int g=0; g<G; g++){ 
      Q[i][g]=0;
      for(unsigned int s=0; s<S; s++){
	Q[i][g]+=pK[s][i][g];
      }
      Q[i][g]/=S;
    }
  }
}


// "Step 2": labelswitch to Q.
int klswitch_step2(double **Q, ERGMM_MCMC_Par *sample, ERGMM_MCMC_Par *tmp, 
		   int S, int n, int d, int G,
		   unsigned int *perm, unsigned int *bestperm, unsigned int *dir,
		   double ***pK){
  int changed=FALSE;
  for(unsigned int s=0; s<S;s++){
    ERGMM_MCMC_Par *cur=sample+s;
    for(unsigned int g=0;g<G;g++){
      perm[g]=g+1;
      dir[g]=0;
    }
    double bestkld=-1;
    do{
      double kld=0;
      for(unsigned int i=0; i<n; i++){
	for(unsigned int g=0; g<G; g++){
	  kld+=pK[s][i][perm[g]-1]*log(pK[s][i][perm[g]-1]/Q[i][g]);
	}
      }
      if(bestkld<-0.5 || kld<bestkld){// Random tiebreaker==bad... || (kld==bestkld && unif_rand()<0.5)){
	if(bestkld>=-0.5) changed=TRUE;
	memcpy(bestperm,perm,G*sizeof(unsigned int));
	bestkld=kld;
      }
    }while(nextperm(G,perm,dir));

    // Now that we have our best permutation...

    apply_perm(bestperm,cur,pK[s],tmp,n,d,G);
    R_CheckUserInterrupt();
  }
  return changed;
}
