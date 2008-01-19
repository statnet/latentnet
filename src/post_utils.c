/********************************************/
/* Utilities for processing posterior draws */
/********************************************/

#include <R.h>
#include <Rmath.h>
#include "post_utils.h"
#include "matrix_utils.h"
#include "ergmm_probs.h"
#include "ergmm_families.h"
#include "mvnorm.h"

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

/* Procustes code due to Raphael Gottardo */

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

R_INLINE void ergmm_par_pred(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par){
  if(model->dir){
    for(unsigned i=0;i<model->verts;i++)
      for(unsigned j=0;j<model->verts;j++)
	model->dY[i][j]+=model->E_edge(model,par,i,j);
  }
  else{
    for(unsigned i=0;i<model->verts;i++)
      for(unsigned j=0;j<=i;j++)
	model->dY[i][j]+=model->E_edge(model,par,i,j);
  }
}

void post_pred_wrapper(int *S, 
		       
		       int *n, int *p, int *d,
		       
		       int *dir,
		       int *family, int *iconsts, double *dconsts,
		       
		       double *vX,
		       
		       double *Z_mcmc, 
		       double *coef_mcmc,
		       double *sender_mcmc, double *receiver_mcmc,
		       
		       int *sociality,
		       int *vobserved_ties,

		       double *vEY,
		       int *s_MKL,
		       int *verbose){
  unsigned int i,j,k;
  unsigned int **observed_ties = (unsigned int **) (vobserved_ties ? Runpack_imatrix(vobserved_ties,*n,*n,NULL) : NULL);
  double ***X = d3array(*p,*n,*n);
  double **Z = dmatrix(*n,*d), *coef = dvector(*p), *sender = sender_mcmc ? dvector(*n):NULL, *receiver = receiver_mcmc ? dvector(*n):NULL;
  
  // set up all of the covariate matrices if covariates are involed 
  // if p=0 (ie no covariates then these next two loops will do nothing)
  //

  for(k=0;k<*p;k++){
    for(i=0;i<*n;i++){
      for(j=0;j<*n;j++){
	X[k][i][j] = vX[ k*(*n)*(*n) + i*(*n) + j ];
      }
    }
  }

  ERGMM_MCMC_Model model = {*dir,
			    NULL, // iY
			    dmatrix(*n,*n), // dY
			    X, // X
			    observed_ties,
			    ERGMM_MCMC_lp_edge[ERGMM_MCMC_to_cont[*family-1]],
			    ERGMM_MCMC_E_edge[ERGMM_MCMC_to_cont[*family-1]],			    
			    0,
			    iconsts,
			    dconsts,
			    *n, // verts
			    *d, // latent
			    *p, // coef
			    0,
			    *sociality};
  
  for(unsigned int s=0; s<*S; s++){
    ERGMM_MCMC_Par par = {*d ? Runpack_dmatrixs(Z_mcmc+s,*n,*d,Z,*S) : NULL, // Z
			  Runpack_dvectors(coef_mcmc+s,*p,coef,*S), // coef
			  NULL, // Z_mean
			  NULL, // Z_var
			  NULL, // Z_pK			  
			  sender_mcmc? Runpack_dvectors(sender_mcmc+s,*n,sender,*S):NULL, // sender
			  0, // sender_var
			  receiver_mcmc? Runpack_dvectors(receiver_mcmc+s,*n,receiver,*S):NULL, // receiver
			  0, // receiver_var
			  NULL, // Z_K
			  0, // llk
			  NULL, // lpedge
			  0, // lpZ		  
			  0, // lpLV
			  0, // lpcoef
			  0 // lpRE
    };
    ergmm_par_pred(&model,&par);
  }

  dscalar_times_matrix(1/(double)*S, model.dY, model.verts, model.verts, model.dY);
  
  Rpack_dmatrixs(model.dY,*n,*n,vEY,1);

  if(s_MKL){
    *s_MKL=-1;
    double m_llk=-HUGE_VAL;
    for(unsigned int s=0; s<*S; s++){
      ERGMM_MCMC_Par par = {*d ? Runpack_dmatrixs(Z_mcmc+s,*n,*d,Z,*S) : NULL, // Z
			    Runpack_dvectors(coef_mcmc+s,*p,coef,*S), // coef
			    NULL, // Z_mean
			    NULL, // Z_var
			    NULL, // Z_pK			  
			    sender_mcmc? Runpack_dvectors(sender_mcmc+s,*n,sender,*S):NULL, // sender
			    0, // sender_var
			    receiver_mcmc? Runpack_dvectors(receiver_mcmc+s,*n,receiver,*S):NULL, // receiver
			    0, // receiver_var
			    NULL, // Z_K
			    0, // llk
			    NULL, // lpedge
			    0, // lpZ		  
			    0, // lpLV
			    0, // lpcoef
			    0 // lpRE
      };
      
      double llk=ERGMM_MCMC_lp_Y(&model,&par,FALSE);

      if(llk>m_llk){
	*s_MKL=s;
	m_llk=llk;
      }
    }
  }

  P_free_all();

  /* Memory freed by GC. */

  return;
}

void klswitch_wrapper(int *maxit, int *S, int *n, int *d, int *G,
		      double *vZ_mcmc, int *Z_ref, double *vZ_mean_mcmc, double *vZ_var_mcmc,
		      int *vZ_K_mcmc, double *vZ_pK_mcmc,
		      double *vQ, int *verbose){

  /* R function enabling uniform RNG (for tie-beaking) */
  GetRNGstate();
  

  if(*verbose>1) Rprintf("KLswitch: Allocating memory.\n");
  ERGMM_MCMC_Par *samples = (ERGMM_MCMC_Par *) P_alloc(*S,sizeof(ERGMM_MCMC_Par));
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
    ERGMM_MCMC_Par *cur=samples+s;
    if(*Z_ref) cur->Z = Z_space[0];
    else cur->Z = Runpack_dmatrixs(vZ_mcmc+s,*n,*d,Z_space[s],*S);
    cur->Z_mean = Runpack_dmatrixs(vZ_mean_mcmc+s,*G,*d,Z_mean_space[s],*S);
    cur->Z_var = Runpack_dvectors(vZ_var_mcmc+s,*G,Z_var_space[s],*S);
    cur->Z_pK = Runpack_dvectors(vZ_pK_mcmc+s,*G,Z_pK_space[s],*S);
    cur->Z_K = Runpack_ivectors(vZ_K_mcmc+s,*n,Z_K_space[s],*S);

    /* Precalculate p(K|Z,mu,sigma).
       This is a speed-memory trade-off:
       we allocate approx S*n*G*8 bytes so that we don't have to do
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
    if(!klswitch_step2(Q, samples, &tmp, *S, *n, *d, *G, perm, bestperm, dir, pK) && it>0){
      if(*verbose>1) Rprintf("KLswitch: Iterating: Q converged after %u iterations.\n", it+1);
      break;
    }
    klswitch_step1(samples, *S, *n, *d, *G, Q, pK);


    if(*verbose>2) Rprintf("KLswitch: Iterating: Completed %u/%d.\n",it+1,*maxit);
  }

  if(*verbose>1) Rprintf("KLswitch: Packing for return to R.\n");

  for(unsigned int s=0; s<*S; s++){
    ERGMM_MCMC_Par *cur=samples+s;
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
R_INLINE int nextperm(unsigned int n, unsigned int *p, unsigned int *dir){
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


R_INLINE void apply_perm(unsigned int *perm, ERGMM_MCMC_Par *to, double **pK, ERGMM_MCMC_Par *tmp, int n, int d, int G){
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
void klswitch_step1(ERGMM_MCMC_Par *samples, int S, int n, int d, int G, double **Q, double ***pK){
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
int klswitch_step2(double **Q, ERGMM_MCMC_Par *samples, ERGMM_MCMC_Par *tmp, 
		   int S, int n, int d, int G,
		   unsigned int *perm, unsigned int *bestperm, unsigned int *dir,
		   double ***pK){
  int changed=FALSE;
  for(unsigned int s=0; s<S;s++){
    ERGMM_MCMC_Par *cur=samples+s;
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
