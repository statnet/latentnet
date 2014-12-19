/*  File src/matrix_utils.h in package latentnet, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2014 Statnet Commons
 */
#ifndef MATRIX_UTILS_H
#define MATRIX_UTILS_H
#include<R.h>
#include"P_alloc.h"

#define FOUND 0
#define NOTFOUND !0

double *dvector(unsigned int n);
int *ivector(unsigned int n);
double **dmatrix(unsigned int n,unsigned int m);
int **imatrix(unsigned int n,unsigned int m);
double ***d3array(unsigned int n1,unsigned int n2, unsigned int n3);

double **Runpack_dmatrix(double *vA, unsigned int n, unsigned int m, double **Aspace);
double **Runpack_dmatrixs(double *vA, unsigned int n, unsigned int m, double **Aspace, unsigned int sample_size);

int **Runpack_imatrix(int *vA, unsigned int n, unsigned int m, int **Aspace);
void Rpack_dmatrixs(double **A, unsigned int n, unsigned int m, double *to, unsigned int sample_size);
void Rpack_d3array(double ***A, unsigned int n1, unsigned int n2, unsigned int n3, double *to);
double *** Runpack_d3array(double *vA, unsigned int n1, unsigned int n2, unsigned int n3, double ***A);
void Rpack_dvectors(double *a, unsigned int n, double *to, unsigned int sample_size);
void Rpack_ivectors(int *a, unsigned int n, int *to, unsigned int sample_size);

void print_dvector (double *a, unsigned int length, FILE *stream );
void print_drvector(double *a, unsigned int n, FILE *stream);
void print_ivector (int *a, unsigned int length, FILE *stream );
void print_dmatrix (double **a, unsigned int nrow, unsigned int ncol, FILE* stream);
void print_imatrix (int **a, unsigned int nrow, unsigned int ncol, FILE* stream);
void init_dvector(double *x, unsigned int n, double value);
void dmatrix_init(double **A, unsigned int n, unsigned int m, double value);
void init_ivector(int *x, unsigned int n, int value);
double *cat_dmectors(double *x, unsigned int nx, double *y, unsigned int ny);
double *cat_dmector_scalar(double *x, unsigned int nx, double y, unsigned int end);
double *dvector_times_matrix(double *x, unsigned int n,double **A, unsigned int m, double *b);
void dscalar_times_matrix(double x, double **A, unsigned int n, unsigned int m, double **B);
void dmatrix_plus_scalar_times_matrix(double x, double **A, unsigned int n, unsigned int m, double **B);
void dmatrix_multiply(double **A,unsigned int na,unsigned int ma, double **B, unsigned int mb, 
		      double **C);
void dmatrix_crossprod(double **A,unsigned int na,unsigned int ma, double **B, unsigned int mb, 
		       double **C);
int *Runpack_ivectors(int *va, unsigned int n, int *a, unsigned int sample_size);
void imatrix_multiply(int **A,unsigned int na,unsigned int ma, int **B, unsigned int mb, int **C);
void dmatrix_addition(double **A, unsigned int n, unsigned int m, double **B);
void imatrix_init(int **A, unsigned int n, unsigned int m, int value);
void t(double **A, unsigned int n, unsigned int m, double **tA);
void dmatrix_copy_contents(double **source,double **dest,unsigned int n,unsigned int m);
double *copy_dvector(double *source,double *dest,unsigned int n);
int *copy_ivector(int *source,int *dest, unsigned int n);
double mean(double *x, unsigned int n);
double *Runpack_dvectors(double *va, unsigned int n, double *a, unsigned int sample_size);

double dmatrix_scale_to(double **A,unsigned int n, unsigned int m, double rms_wanted);
void dmatrix_scale_by(double **A, unsigned int n, unsigned int m, double by);
void dvector_scale_by(double *v, unsigned int n, double by);
int dgesvd_full_wrapper(double **A, int n, int m, double **U, double *S, double **tV, double *workspace, int lworkspace);
#endif
