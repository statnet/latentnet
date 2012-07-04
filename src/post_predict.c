/**********************************************************/
/* Utilities for computing posterior dyad expected values */
/**********************************************************/

#include <R.h>
#include <Rmath.h>
#include "post_predict.h"
#include "matrix_utils.h"
#include "ergmm_probs.h"
#include "ergmm_families.h"
#include "ergmm_latent_effects.h"

/*R_INLINE*/ void ergmm_par_pred(ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par){
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
		       
		       int *n, int *p, int *d, int *latent_eff, int *family, int *res,
		       
		       int *dir,
		       int *iconsts, double *dconsts,

		       
		       double *vX,
		       
		       double *Z_mcmc, 
		       double *coef_mcmc,
		       double *sender_mcmc, double *receiver_mcmc,
		       
		       int *vobserved_ties,

		       double *vEY,
		       int *s_MKL,
		       int *verbose){

  // This was added because newer versions of R no longer pass a 0-length vector as NULL, so we have to do it here.
  if(*p==0){
    vX=coef_mcmc=NULL;
  }
  if(*d==0){
    Z_mcmc=NULL;
  }
  if(res[0]==0&&res[2]==0){
    sender_mcmc=NULL;
  }
  if(res[1]==0){
    receiver_mcmc=NULL;
  }


  unsigned int i,j,k;
  unsigned int **observed_ties = (unsigned int **) (*vobserved_ties>=0 ? Runpack_imatrix(vobserved_ties,*n,*n,NULL) : NULL);
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
			    res[2],
			    latent_eff ? ERGMM_MCMC_latent_eff[*latent_eff-1] : NULL
  };
  
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
    if(model.sociality) par.receiver=par.sender;
    ergmm_par_pred(&model,&par);
  }

  dscalar_times_matrix(1/(double)*S, model.dY, model.verts, model.verts, model.dY);
  
  Rpack_dmatrixs(model.dY,*n,*n,vEY,1);

  if(*s_MKL){
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

