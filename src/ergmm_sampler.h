#ifndef ERGMM_SAMPLER_H
#define ERGMM_SAMPLER_H

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
   [1] RE_delta
   [2] Z_scl_delta
   [3]- group_deltas
*/
#define GROUP_DELTAS_START 3

void ERGMM_MCMC_wrapper(int *samples_stored, int *interval,
			   
			int *n, int *p, int *d, int *G,
			  
			int *dir, int *viY, double *vdY,
			int *family, int *iconsts, double *dconsts,

			double *vX,
			  
			double *llk_mcmc, double *lpZ_mcmc, double *lpcoef_mcmc, double *lpRE_mcmc, double *lpLV_mcmc, double *lpREV_mcmc,
			   
			double *vZ_start,

			double *Z_pK_start, double *vZ_mean_start,double *Z_var_start,int *Z_K_start,

			double *Z_var_prior, double *Z_mean_prior_var, 
			double *Z_K_prior, double *Z_var_prior_df,

			double *Z_mcmc, double *Z_rate_move, double *Z_rate_move_all,

			int *Z_K_mcmc, double *Z_pK_mcmc, double *Z_mean_mcmc, double *Z_var_mcmc,
			  
			double *coef_start,
			double *coef_mean, double *coef_prior_var,
			double *coef_mcmc, double *coef_rate, 
			  
			  
			double *sender_start, double *receiver_start,
			double *sender_var_start, double *receiver_var_start,

			double *sender_var_prior, double *sender_var_prior_df,
			double *receiver_var_prior, double *receiver_var_prior_df,

			double *sender_mcmc, double *receiver_mcmc,
			double *sender_var_mcmc, double *receiver_var_mcmc,

			int *sociality,
			int *vobserved_ties,
			double *deltas);

void ERGMM_MCMC_init(unsigned int samples_stored, unsigned int interval, 

		     unsigned int n, 
		     unsigned int p, unsigned int d, unsigned int G,

		     unsigned int dir, int **iY, double **dY,

		     unsigned int family, int *iconsts, double *dconsts,

		     double ***X,

		     double *llk_mcmc, double *lpZ_mcmc, double *lpcoef_mcmc, double *lpRE_mcmc, double *lpLV_mcmc, double *lpREV_mcmc,

		     double **Z_start,
		     double *Z_pK_start, double **Z_mean_start, double *Z_var_start, unsigned int *Z_K_start,
		     double Z_var_prior, double Z_mean_prior_var, double Z_K_prior,
		     double Z_var_prior_df,
		     double *Z_mcmc, double *Z_rate_move, double *Z_rate_move_all, int *K_mcmc,
		     double *Z_pK_mcmc,
		     double *Z_mean_mcmc, double *Z_var_mcmc,

		     double *coef_start,
		     double *coef_mcmc, double *coef_rate, 
		     double *coef_mean, double *coef_prior_var, 

		     double *sender_start, double *receiver_start,
		     double sender_var_start, double receiver_var_start,
		     double sender_var_prior, double sender_var_prior_df,
		     double receiver_var_prior, double receiver_var_prior_df,
		     double *sender_mcmc, double *receiver_mcmc,
		     double *sender_var_mcmc, double *receiver_var_mcmc,
		     unsigned int sociality,
		     unsigned int **observed_ties,

		     double Z_delta,
		     double RE_delta,
		     double Z_scl_delta,
		     double **group_deltas, unsigned int group_prop_size);

void ERGMM_MCMC_loop(ERGMM_MCMC_Model *model, ERGMM_MCMC_Priors *prior,
		     ERGMM_MCMC_MCMCState *cur, ERGMM_MCMC_MCMCSettings *setting,
		     ERGMM_MCMC_ROutput *outlists);

void ERGMM_MCMC_store_iteration(unsigned int pos, ERGMM_MCMC_Model *model, ERGMM_MCMC_Par *par,
		     ERGMM_MCMC_MCMCSettings *setting, ERGMM_MCMC_ROutput *outlists);
#endif
