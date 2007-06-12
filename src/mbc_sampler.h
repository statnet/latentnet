#ifndef MBC_SAMPLER_H
#define MBC_SAMPLER_H 1

#include "ergmm_structs.h"
#define FALSE 0
#define TRUE !0

/* First 1 positions in the outlists are reserved for special values:
   [0] Iteration with the highest likelihood so far.
*/
#define MBC_OUTLISTS_RESERVE 1
void MBC_MCMC_wrapper(int *samples_stored,
		      int *interval,

		      int *n,
		      int *d,
		      int *G,
			  
		      double *lpZList, 
		      double *lpLVList, 
		      
		      double *vZ,
		      
		      double *epsilon,
		      double *mu,
		      double *Sigma,
		      int *Ki,
		      
		      double *Sigprior, 
		      double *muSigprior, 
		      double *dirprior,
		      double *alphaprior,
		      
		      int *KiList, 
		      double *Z_pKList, 
		      double *muList, 
		      double *SigmaList);
void MBC_MCMC_init(unsigned int samples_stored, 
		   unsigned int interval, 

		   unsigned int n,
		   unsigned int d,
		   unsigned int G,
		   
		   double *lpZList,
		   double *lpLVList, 
		   
		   double **Z,

		   double *epsilon, 
		   double **Z_mu_start, 
		   double *Sigma, 
		   unsigned int *Ki,

		   double Sigprior,
		   double muSigprior,
		   double dirprior,
		   double alphaprior,

		   int *KList,
		   double *Z_pKList,
		   double *muList,
		   double *SigmaList);

void MBC_MCMC_loop(ERGMM_MCMC_Model *model, ERGMM_MCMC_Priors *prior,
		   ERGMM_MCMC_MCMCState *cur, ERGMM_MCMC_MCMCSettings *setting, ERGMM_MCMC_ROutput *outlists);

void MBC_MCMC_store_iteration(unsigned int pos, ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par,
			      ERGMM_MCMC_MCMCSettings *setting, ERGMM_MCMC_ROutput *outlists);

#endif
