/****************************************/
/* Sampling from Dirichlet distribution */
/****************************************/

#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "matrix_utils.h"


/* Dirichlet generator */
void rdirich(unsigned int n, double *epsilon) 
{
  unsigned int i;
  double  sum = 0;
  //  Rprintf("n = %d\n",n);
  for (i=0; i<n; ++i){
    //    if(epsilon[i]==1)epsilon[i] +=0.5;
    //    Rprintf("epsilon[%d] = %1.7f\n",i,epsilon[i]);
    epsilon[i] = rgamma(epsilon[i], 1.0);
    sum += epsilon[i];
  }
  
  for (i=0; i<n; ++i)
    epsilon[i] = epsilon[i]/sum;
}

void dirichlet_wrapper(int *n, double *epsilon)
{
  GetRNGstate(); 
  rdirich(*n, epsilon);
  PutRNGstate();
}

