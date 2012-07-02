#ifndef POST_PREDICT_H
#define POST_PREDICT_H

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
		       int *verbose);
#endif
