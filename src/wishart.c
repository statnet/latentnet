/*  File src/wishart.c in package latentnet, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2018 Statnet Commons
 */
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

