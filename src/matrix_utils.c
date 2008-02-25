/****************************************/
/* Matrix, vector, and memory utilities */
/****************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include "matrix_utils.h"


/*  Allocates memory for a vector of doubles of length n */
double *dvector(unsigned int n){
  if(n<=0) return NULL;

  double *a;
  unsigned int i;
  a = (double*) P_alloc(n,sizeof(double));
  if(a == NULL){
    P_free_all();
    error("Not enough memory to make double vector.");
  }

  for(i=0;i<n;i++){
    a[i]=0.0;
  }
  return a;
}

/*  Allocates memory for a vector of doubles of length n */
int *ivector(unsigned int n){
  if(n<=0) return NULL;

  int *a;
  unsigned int i;
  a = (int*) P_alloc(n,sizeof(int));
  if(a == NULL){
    P_free_all();
    error("Not enough memory to make integer vector.");
  }

  for(i=0;i<n;i++){
    a[i]=0;
  }
  return a;
}

/*  Allocates memory for an n by m matrix of doubles */
double **dmatrix(unsigned int n,unsigned int m)
{
  if(n<=0 || m<=0) return NULL;

  double **A;
  unsigned int i, j;

  /* assigning memory and initialize */
  A = (double**) P_alloc(n,sizeof(double*));
  if(A == NULL){
    P_free_all();
    error("Not enough memory to make double matrix.");
  }
  A[0] = (double *) P_alloc(n*m,sizeof(double));
  if(A[0] == NULL){
    P_free_all();
    error("Not enough memory to make double matrix.");
  }
  for(i=0;i<n;i++){
    A[i] = A[0]+i*m;
    for(j=0;j<m;j++){
      A[i][j]=0.0;
    }
  }

  return A;
}

/*  Allocates memory for an n1 by n2 by n3 matrix of doubles */
double ***d3array(unsigned int n1,unsigned int n2, unsigned int n3)
{
  if(n1<=0 || n2<=0 || n3<=0) return NULL;

  double ***A;

  /* allocate space for the matrix and the pointers */
  A = (double***) P_alloc(n1,sizeof(double**));
  if(A == NULL){
    P_free_all();
    error("Not enough memory to make 3D array.");
  }
  A[0] = (double **) P_alloc(n1*n2,sizeof(double*));
  if(A[0] == NULL){
    P_free_all();
    error("Not enough memory to make 3D array.");
  }
  A[0][0] = (double *) P_alloc(n1*n2*n3,sizeof(double));
  if(A[0] == NULL){
    P_free_all();
    error("Not enough memory to make 3D array.");
  }

  
  for(unsigned int i1=0;i1<n1;i1++){
    A[i1] = A[0]+i1*n2;
    A[i1][0] = A[0][0]+i1*n2*n3;
    for(unsigned int i2=1;i2<n2;i2++){
      A[i1][i2] = A[i1][0]+i2*n3;
      for(unsigned int i3=0;i3<n3;i3++){
	A[i1][i2][i3]=0.0;
      }
    }
  }

  return A;
}

/*  Allocates memory for an n by m matrix of doubles */
int **imatrix(unsigned int n,unsigned int m)
{
  if(n<=0 || m<=0) return NULL;

  int **A;
  unsigned int i, j;
  /* assigning memory and initializing */
  A = (int**) P_alloc(n,sizeof(int*));
  if(A == NULL){
    P_free_all();
    error("Not enough memory to make integer matrix.");
  }
  A[0] = (int*) P_alloc(n*m,sizeof(int));
  if(A[0] == NULL){
    P_free_all();
    error("Not enough memory to make integer matrix.");
  }
  for(i=0;i<n;i++){
    A[i] = A[0]+i*m;
    for(j=0;j<m;j++){
      A[i][j]=0;
    }
  }
  return A;
}

/* Deserializes an array from R into a matrix.
   Note that the "top" offset needs to already be applied to "to".
 */
double **Runpack_dmatrix(double *vA, unsigned int n, unsigned int m, double **Aspace){
  if(!Aspace) Aspace = dmatrix(n,m); 
  unsigned int i,j;
  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      Aspace[i][j]=vA[i+j*n];
  return(Aspace);
}

double **Runpack_dmatrixs(double *vA, unsigned int n, unsigned int m, double **Aspace, unsigned int sample_size){
  if(!Aspace) Aspace = dmatrix(n,m); 
  unsigned int i,j;
  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      Aspace[i][j] = vA[(i+n*j)*sample_size];
  return Aspace;
}
/* Deserializes a matrix from R into a vector.
   Note that the "top" offset needs to already be applied to "va".
 */
double *Runpack_dvectors(double *va, unsigned int n, double *a, unsigned int sample_size){
  if(!a) a=dvector(n);
  unsigned int i;
  for(i=0;i<n;i++)
    a[i]=va[sample_size*i];
  return a;
}

int *Runpack_ivectors(int *va, unsigned int n, int *a, unsigned int sample_size){
  if(!a) a=ivector(n);
  unsigned int i;
  for(i=0;i<n;i++)
    a[i]=va[sample_size*i];
  return a;
}

/* Deserializes an (integer) array from R into a matrix.
   Note that the "top" offset needs to already be applied to "to".
 */
int **Runpack_imatrix(int *vA, unsigned int n, unsigned int m, int **Aspace){
  if(!Aspace) Aspace = imatrix(n,m); 
  unsigned int i,j;
  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      Aspace[i][j]=vA[i+j*n];
  return(Aspace);
}

/* Serializes a matrix for returning to R.
   Note that the "top" offset needs to already be applied to "to".
 */
void Rpack_dmatrixs(double **A, unsigned int n, unsigned int m, double *to, unsigned int sample_size){
  unsigned int i,j;
  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      to[(i+n*j)*sample_size] = A[i][j];
}

/* Serializes a dvector into an array for returning to R.
   Note that the "top" offset needs to already be applied to "to".
 */
void Rpack_dvectors(double *a, unsigned int n, double *to, unsigned int sample_size){
  unsigned int i;
  for(i=0;i<n;i++)
    to[sample_size*i]=a[i];
}

/* Serializes an ivector into an array for returning to R.
   Note that the "top" offset needs to already be applied to "to".
 */
void Rpack_ivectors(int *a, unsigned int n, int *to, unsigned int sample_size){
  unsigned int i;
  for(i=0;i<n;i++)
    to[sample_size*i]=a[i];
}

// free_* functions were here. I have switched to R's garbage collection, so they are not needed.

/* prints to file stream*/
void print_dvector(double *a, unsigned int n, FILE *stream) {
  unsigned int i;
  for (i=0;i<n;i++) {
    // printf("hidie\n");
    fprintf(stream, "%+.8lf ", a[i]);
    //fprintf(stream, "%+.8e ", a[i]);
  }
  fprintf(stream, "\n");
}

void print_drvector(double *a, unsigned int n, FILE *stream) {
  unsigned int i;
  for (i=0;i<n;i++) {
    fprintf(stream, "%e\n", a[i]);
  }
  fprintf(stream, "\n");
}

/* prints to file stream*/
void print_ivector(int *a, unsigned int n, FILE *stream) {
  unsigned int i;
  for (i=0;i<n;i++) {
    if(a[i]<10)
      fprintf(stream, "%d   ", a[i]);
    else if(a[i]<100)
      fprintf(stream, "%d  ", a[i]);
//    else
  //    fprintf(stream, "%d  ", a[i]);
  }
  fprintf(stream, "\n");
}

/* prints to file stream*/
void print_dmatrix(double **a, unsigned int n, unsigned int m, FILE* stream) {
  unsigned int i;
  //printf("hello in print_dmatrix.  n is: %d, m is %d\n",n,m);  
  //printf("here is the address of a: %e",a);
  for (i=0;i<n;i++) {
    //   printf("hi\n");
    print_dvector (a[i], m, stream);
  }
}

/* prints to file stream*/
void print_imatrix (int **a, unsigned int n, unsigned int m, FILE* stream) {
  unsigned int i;
  for (i=0;i<n;i++) {
    print_ivector (a[i], m, stream);
  }
}

void init_dvector(double *x, unsigned int n, double value)
{
  unsigned int i;
  for(i=0;i<n;i++){
    x[i] = value;
  }
}

void init_dmatrix(double **A, unsigned int n, unsigned int m, double value)
{
  unsigned int i,j;
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      A[i][j]=value;
    }
  }
}

void init_ivector(int *x, unsigned int n, int value)
{
  unsigned int i;
  for(i=0;i<n;i++){
    x[i] = value;
  }
}



/* Concatenates two vectors, the second on the bottom */
double *cat_dvectors(double *x, unsigned int nx, double *y, unsigned int ny)
{
  unsigned int i;
  double *xy = dvector(nx+ny);

  for(i=0;i<nx;i++){
    xy[i] = x[i];
  }
  for(i=0;i<ny;i++){
    xy[i+nx] = y[i];
  }
  return(xy);
}

/* Concatenates a vector and a scalar */
/* if end equals true then place at end */
double *cat_dvector_scalar(double *x, unsigned int nx, double y, unsigned int end)
{
  unsigned int i;
  double *xy; 
  xy = dvector(nx+1);

  if(end==0){
    for(i=0;i<nx;i++){
      xy[i] = x[i];
    }
    xy[nx] = y;
  }
  else{
    xy[0] = y;
    for(i=1;i<=nx;i++){
      xy[i] = x[i-1];
    }
  }
  return(xy);
}


/* Computes AX, where A is na by ma and B is nb by mb */
/*  b needs ot be a vector of length m*/
double *dvector_times_matrix(double *x, unsigned int n, double **A, unsigned int m, double *b)
{
  
  unsigned int i=0,j=0;
  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      b[i] += x[j]*A[j][i];
    }
  }
  return(b);
}

/* mutliply the scalar times the matrix and return the result in B*/
void dscalar_times_matrix(double x, double **A, unsigned int n, unsigned int m, double **B)
{
  unsigned int i,j;
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      B[i][j] = A[i][j] * x;
    } 
  }
}

/* multiply the scalar times the matrix and add the result to B*/
void dmatrix_plus_scalar_times_matrix(double x, double **A, unsigned int n, unsigned int m, double **B)
{
  unsigned int i,j;
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      B[i][j] += A[i][j] * x;
    } 
  }
}

/* Computes C0+AB, where A is na by ma and B is ma by mb, and C0 is the initial contents of C. */
void dmatrix_multiply(double **A,unsigned int na,unsigned int ma, double **B, unsigned int mb, 
			  double **C)
{
  unsigned int i=0,j=0,k=0;
 
  for(i=0;i<na;i++){
    for(j=0;j<mb;j++){
      for(k=0;k<ma;k++){
	C[i][j] += A[i][k]*B[k][j];
      }
    }
  }
  return;
}

/* Computes AB, where A is na by ma and B is ma by mb */
void imatrix_multiply(int **A,unsigned int na,unsigned int ma, int **B, unsigned int mb, int **C)
{
  unsigned int i=0,j=0,k=0;
 
  /* if(ma != nb) return (double**)(-1); */
  for(i=0;i<na;i++){
    for(j=0;j<mb;j++){
      for(k=0;k<ma;k++){
	C[i][j] += A[i][k]*B[k][j];
      }
    }
  }
  return;
}

/* add the second matrix to the first */
void dmatrix_addition(double **A, unsigned int n, unsigned int m, double **B)
{
  unsigned int i,j;
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      A[i][j] += B[i][j];
    }
  }
}


/* Makes sure that all values of the matrix are the same*/
void init_imatrix(int **A, unsigned int n, unsigned int m, int value)
{
  unsigned int i,j;

  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      A[i][j] = value;
  return;
}


/* Returns the transpose of the matrix */
void t(double **A, unsigned int n, unsigned int m, double **tA)
{
   unsigned int i,j;
   for(i=0;i<n;i++){
     for(j=0;j<m;j++){
       tA[j][i] = A[i][j];
     }  
   }
   return;
}


void copy_dmatrix(double **source,double **dest,unsigned int n,unsigned int m)
{
  unsigned int i;
  for(i=0;i<n;i++)
    memcpy(dest[i],source[i],m*sizeof(double));
  return;
}

double* copy_dvector(double *source,double *dest,unsigned int n){
  if(!dest) dest=dvector(n);
  memcpy(dest,source,n*sizeof(double));
  return dest;
}

int* copy_ivector(int *source,int *dest,unsigned int n){
  if(!dest) dest=ivector(n);
  memcpy(dest,source,n*sizeof(int));
  return dest;
}

double mean(double *x, unsigned int n){
  unsigned int i=0;
  double mu=0.0;
  for(i=0;i<n;i++){
    mu += x[i];
  }
  return(mu/n);
}

/*R_INLINE*/ double dvector_dist(double *u, double *v, unsigned int dim){
  unsigned int k;
  double dist,dist2=0;
  for(k=0;k<dim;k++){
    dist=u[k]-v[k];
    dist2+=dist*dist;
  }
  return(sqrt(dist2));
}

double dmatrix_scale_to(double **A,unsigned int n, unsigned int m, double rms_wanted){
  double rms=0;
  unsigned int i,j;
  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      rms+=A[i][j]*A[i][j];
  rms=sqrt(rms/(n*m));
  dmatrix_scale_by(A,n,m,rms_wanted/rms);
  return(rms_wanted/rms);
}

void dmatrix_scale_by(double **A, unsigned int n, unsigned int m, double by){
  unsigned int i,j;
  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      A[i][j]*=by;
}

void dvector_scale_by(double *v, unsigned int n, double by){
  unsigned int i;
  for(i=0;i<n;i++)
    v[i]*=by;
}

/* Inverts a matrix "in place".
   "workspace" must be an expendable block of memory of length (n*n*3+4*n)*sizeof(double)*/
/*R_INLINE*/ int inverse(double **x, int n, double **res, double *workspace)
{
  /** QR decomposition solve x*x=I **/
  /** res contains the result **/

  int i,j, info = 0, rank, *pivot, p;
  double tol = 1.0E-7, *qraux, *work;
  double *xt, *yt, *rest;
  
  qraux = workspace;
  work  = qraux+n;
  rest = work+2*n;
  yt = rest+n*n;
  xt = yt+n*n;
  pivot = (int *) (xt+n*n);
 
  for(i = 0; i < n; i++)
    pivot[i] = i+1;
  
  /** Copy the matrix by column **/
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      xt[i*(n)+j]=x[j][i];
  

  p = n;
  
  F77_CALL(dqrdc2)(xt, &n, &n, &p, &tol, &rank,qraux, pivot, work);
  
  
  /** Copy the matrix by column **/

  for(i=0;i<n;i++)
      for(j=0;j<n;j++)
	yt[i*(n)+j]= j==i?1:0;
  
  F77_CALL(dqrcf)(xt, &n, &rank, qraux,yt, &n, rest, &info);
  
  /** Put back into a matrix **/
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      res[j][i]=rest[i*(n)+j];
  return(rank==n? FOUND:NOTFOUND);
}

/* vectors is non zero if want eigen vectors returned.
   length(EValues) = n
   dim(EVectors) = n,n
   "workspace" must be an expendable block of memory of length (n*n*2+n*3)*sizeof(double).
*/
/*R_INLINE*/ int sym_eigen(double **A, int n, int vectorsflag, double *EValues, double **EVectors,double *workspace)
{
  int err=0,i=0,j=0;
  /*
  double *vA = dvector(n*n);
  double *vEVectors = dvector(n*n);
  double *l=dvector(n);
  double *fv1 = dvector(n);
  double *fv2 = dvector(n);
  */
  //  Workaround...
  double *vA = workspace;
  double *vEVectors = vA+n*n;
  double *l=vEVectors+n*n;
  double *fv1 = l+n;
  double *fv2 = fv1+n;

  /** Copy the matrix by column **/
  /* make A into a vector to pass in*/
  for(j=0;j<n;j++){
    for(i=0;i<n;i++){
      vA[(n*j)+i] = A[i][j];
    }
  }
  
  
  /*  int F77_NAME(rs)(int *nm, int *n, double *a, double *w,
      int *matz, double *z, double *fv1, double *fv2, int *ierr) */
  F77_NAME(rs)(&n, &n, vA, l, &vectorsflag, vEVectors, fv1, fv2, &err);
  
  for(i=0;i<n;i++){
    EValues[i] = l[n-1-i];
  }  
  
  /* put A and EVectors into matrices to return */
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      EVectors[i][j] = vEVectors[i + (n-(j+1))*n];
    }
  } 
  
  return 0;
}

