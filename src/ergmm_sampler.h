#ifndef ERGMM_SAMPLER_H
#define ERGMM_SAMPLER_H 1

/****************************************************************************/
/*  Original Author: Susan Shortreed, susanms@stat.washington.edu           */
/*  Updated by: Jeremy Tantrum, tantrum@stat.washington.edu                 */
/*  Rewritten by: Pavel Krivitsky, pavel@stat.washington.edu                */
/*  Purpose: main functions for parameter estimation model 2                */
/*           All of this code is for an R function which is incorporated    */
/*           into the R ERGM package.                                       */
/****************************************************************************/


#include "ergmm_structs.h"
#define FALSE 0
#define TRUE !0

/* First 1 positions in the outlists are reserved for special values:
   [0] Iteration with the highest likelihood so far.
   [1] Iteration with the highest joint density of all variables (except K) so far.
*/
#define ERGMM_OUTLISTS_RESERVE 2
#define GET_DEFAULT(p,i,d) ((p)?(p)[(i)]:(d))


/* deltas have the following values:
   [0] Z_delta
   [1] Z_tr_delta
   [2] Z_scl_delta 
   [3] RE_delta
   [4] RE_shift_delta
   [5]-[5+p-1] coef_delta
*/
#define COEF_DELTA_START 5

void ERGMM_MCMC_wrapper(int *samples_stored, int *interval,
			   
			int *n, int *p, int *d, int *G,
			  
			int *dir, int *vY,
			int *family, int *iconsts, double *dconsts,

			double *vX,
			  
			double *llkList, double *lpZList, double *lpcoefList, double *lpREList, double *lpLVList, double *lpREVList,
			   
			double *vZ_start,

			double *epsilon, double *mu,double *Sigma,int *Ki,

			double *Sigprior, double *muSigprior, 
			double *dirprior, double *alphaprior,

			double *vZ_post, double *Z_rate_move, double *Z_rate_move_all,

			int *KiList, double *Z_pKList, double *muList, double *SigmaList,
			  
			double *coef_start,
			double *coef_mean, double *coef_var,
			double *Coef, double *B_rate, 
			  
			  
			double *sender, double *receiver,
			double *sender_var, double *receiver_var,

			double *sender_var_prior, double *sender_var_prior_df,
			double *receiver_var_prior, double *receiver_var_prior_df,

			double *senderList, double *receiverList,
			double *sender_varList, double *receiver_varList,

			int *sociality,
			int *vobserved_ties,
			double *deltas);

void ERGMM_MCMC_init(unsigned int samples_stored, unsigned int interval, 

		     unsigned int n, 
		     unsigned int p, unsigned int d, unsigned int G,

		     unsigned int dir, int **Y,

		     unsigned int family, int *iconsts, double *dconsts,

		     double ***X,

		     double *llkList, double *lpZList, double *lpcoefList, double *lpREList, double *lpLVList, double *lpREVList,

		     double **Z_start,
		     double *epsilon, double **Z_mu_start, double *Sigma, unsigned int *Ki,
		     double Sigprior, double muSigprior, double dirprior,
		     double alphaprior,
		     double *ZList, double *Z_rate_move, double *Z_rate_move_all, int *KList,
		     double *Z_pKList,
		     double *muList, double *SigmaList,

		     double *coef_mle,
		     double *coefList, double *coef_rate, 
		     double *coef_mean, double *coef_var, 

		     double *sender, double *receiver,
		     double sender_var, double receiver_var,
		     double sender_var_prior, double sender_var_prior_df,
		     double receiver_var_prior, double receiver_var_prior_df,
		     double *senderList, double *receiverList,
		     double *sender_varList, double *receiver_varList,
		     unsigned int sociality,
		     unsigned int **observed_ties,

		     double Z_delta, double Z_tr_delta, double Z_scl_delta,
		     double RE_delta, double RE_shift_delta,
		     double *coef_delta);

void ERGMM_MCMC_loop(ERGMM_MCMC_Model *model, ERGMM_MCMC_Priors *prior,
		     ERGMM_MCMC_MCMCState *cur, ERGMM_MCMC_MCMCSettings *setting,
		     ERGMM_MCMC_ROutput *outlists);

void ERGMM_MCMC_store_iteration(unsigned int pos, ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par,
		     ERGMM_MCMC_MCMCSettings *setting, ERGMM_MCMC_ROutput *outlists);
#endif
