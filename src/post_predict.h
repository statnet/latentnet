#ifndef POST_PREDICT_H
#define POST_PREDICT_H

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
		       int *verbose);
#endif
