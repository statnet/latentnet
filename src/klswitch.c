/****************************************************************************/
/*  Author: Jeremy Tantrum, tantrum@stat.washington.edu                     */
/*  Purpose: More support functions for model 2                             */
/*           proposed by Adrian E. Raftery, Mark S. Handcock and Jeremy T   */
/****************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>

//#include "latentUtil.h"

#define ARRAY(x1,x2,x3,n1,n2,n3) ((x1) + ((n1)*(x2)) + ((n1)*(n2)*(x3)))
#define MAT(row,col,numrow) ((row)+((numrow)*(col)))

void klswitch(long p, long nnodes, long nsamp, long ngroups, 
	      long npermute, double *Z, double *mu, double *Sigma, 
	      double *qig, long *permute, long *minat);

void call_klswitch(int *p, int *nnodes, int *nsamp, int *ngroups, 
		   int *npermute, double *Z, double *mu, double *Sigma, 
		   double *qig, int *permute, int *minat)
{
  klswitch((long)*p, (long)*nnodes, (long)*nsamp, (long)*ngroups,
	   (long)*npermute, Z, mu, Sigma, qig, 
	   (long*)permute, (long*)minat);
  /*  Note:  The above variables are coerced to long but they start as    */
  /*         int:  It's impossible to pass longs from R to call_klswitch. */
  /*         This could eventually be changed, since there is no benefit  */
  /*         of coercing int to long.                                     */
}

void klswitch(long p, long nnodes, long nsamp, long ngroups, 
	      long npermute, double *Z, double *mu, double *Sigma, 
	      double *qig, long *permute, long *minat)
{
  /* p = 2, nnodes = 69, nsamp = 995, ngroups = 5, npermute = 5! */
  /* Take in: Z[69,2,995], mu[995,2*5], permute[5!,5], Sigma[995,5], qig[69,5] */
  /* Output: minat[995], minat[i] \in 1..ngroups */
  /*   -0.5[ (x-mu)^t Sig^-1 (x-mu) + p * log(2 * Pi * det(Sig)) ] */
  int i,j,k,g, loop;
  double jexp1, jexp2, prob, probmin, probtemp;

  /*   Rprintf("p=%d nnodes=%d nsamp=%d ngroups=%d npermute=%d\n",p,nnodes,nsamp,ngroups,npermute); */
  /*   Rprintf("Z[1,,1] = [%1.4f, %1.4f]\n",Z[ARRAY(0,0,0,nnodes,p,nsamp)],Z[ARRAY(0,1,0,nnodes,p,nsamp)]); */
  /*   Rprintf("mu[1,,1] = [%1.4f, %1.4f]\n",mu[ARRAY(permute[MAT(0,0,npermute)], 0 ,0,ngroups,p,nsamp)],mu[ARRAY(permute[MAT(0,0,npermute)], 0 ,0,ngroups,p,nsamp)]); */

  for(k=0; k< nsamp; k++)
  {
    probmin = 10000000;
    for(loop=0; loop < npermute; loop ++)
    {
      prob = 0.0;
      for(j=0; j < nnodes; j++)
	for(g=0; g < ngroups; g++)
	{
	  jexp1 = 0.0;
	  for(i=0; i < p; i++)
	  {
	    jexp2 = Z[ARRAY(j,i,k,nnodes,p,nsamp)] - 
	      mu[ARRAY(permute[MAT(loop,g,npermute)]-1, i ,k,ngroups,p,nsamp)];
	    jexp1 += jexp2 * jexp2;
	  }
	  jexp1 = jexp1 / (-2 * Sigma[MAT(k,permute[MAT(loop,g,npermute)]-1,nsamp)]);
	  probtemp = -0.5*(p) * 
	   log(2 * M_PI * Sigma[MAT(k,permute[MAT(loop,g,npermute)]-1,nsamp)]) + jexp1;
	  prob = prob + exp(probtemp) * log(exp(probtemp) / qig[MAT(j,g,nnodes)]);
	}/* end for gh in 0:ngroups */
      if(prob < probmin)
      {
	probmin = prob;
	minat[k] = loop;
      }
    }/* end for loop in 0:npermute */
  }/* end for k in 0:nsamp */
  
}







