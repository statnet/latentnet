#ifndef MVNORM_H
#define MVNORM_H

double dindnormmu(unsigned int n, double *x, double *mu, double sigma, int give_log);
double diidnorm0(unsigned int n, double *x, double sigma, int give_log);
double call_mvdlnorm2(int p, double *mu, double *Sigma, double *x);
void mvdlnorm2(int p, double *mu, double *Sigma, double *x, double *result);

double sdlnorm(unsigned int p, unsigned int ng, unsigned int grp, double **mu, double *Sigma, double *x);
#endif
