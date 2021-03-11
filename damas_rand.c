/*==========================================================
 * damas.c - MEX implementation of DAMAS
 *
 *
 * p = damas(A, b, niter)
 *
 *
 *========================================================*/

#include "mex.h"
#include <stdlib.h>

/* The computational routine */
void damas(double* A, double* y, mwSize N, int niter, double* z)
{
    mwSize i;
    mwSize j;


    int iter = 0;
    int k = 0;

    double new = 0;

    for (iter = 0 ; iter < niter; iter++)
    {
      k = (k + (rand() % N))%N;

        new = y[k];
        for (j = 0 ; j < N ; j++)
        {
          new = new - A[j + k*N] * z[j];

        }
          new = new + A[k + k*N] * z[k];

	       z[k] = (new > 0 ?  new : 0) / A[k + k*N];

    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *A;               /* 1xN input matrix */
    double *y;
    double niter;              /* input scalar */

    double *outMatrix;              /* output matrix */

    size_t N;

    N = mxGetN(prhs[0]);

    plhs[0] = mxCreateDoubleMatrix((mwSize)N,1,mxREAL);
    outMatrix = mxGetPr(plhs[0]);


    /* create a pointer to the real data in the input matrix  */
    A = mxGetPr(prhs[0]);
    y = mxGetPr(prhs[1]);
    niter = mxGetScalar(prhs[2]);
    double* x0 = mxGetPr(prhs[3]);

    for (int i=0; i<N; i++) {
      outMatrix[i] = x0[i];
    }


    damas(A,y,N,niter, outMatrix);
}
