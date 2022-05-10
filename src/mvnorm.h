/*  File src/mvnorm.h in package latentnet, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2022 Statnet Commons
 */
#ifndef MVNORM_H
#define MVNORM_H

double dindnormmu(unsigned int n, double *x, double *mu, double sigma, int give_log);
double diidnorm0(unsigned int n, double *x, double sigma, int give_log);
double call_mvdlnorm2(int p, double *mu, double *Sigma, double *x);
void mvdlnorm2(int p, double *mu, double *Sigma, double *x, double *result);

double sdlnorm(unsigned int p, unsigned int ng, unsigned int grp, double **mu, double *Sigma, double *x);
#endif
