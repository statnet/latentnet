/*  File src/post_predict.h in package latentnet, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2022 Statnet Commons
 */
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
           double *dispersion_mcmc,
		       int *vobserved_ties,
		       double *vEY,
		       int *s_MKL,
		       int *verbose);
#endif
