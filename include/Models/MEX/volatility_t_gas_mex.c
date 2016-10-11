#include "mex.h"
#include "math.h"

void volatility_t_gas_mex(double *y, mwSignedIndex N, mwSignedIndex T,
        double *theta, double *h)
{
//    double *omega; 
    mwSignedIndex i, j;
    double *A, *nu_con;
    double tmp;
        
    /* Variable size arrays */
    A = mxMalloc((N)*sizeof(double));
    nu_con = mxMalloc((N)*sizeof(double));
    
    
    /* Initialise */
    for (i=0; i<N; i++)
    {
        h[i] = theta[i+N]/(1-theta[i+3*N]);
        A[i] = theta[i+2*N]*(theta[i+4*N]+3)/theta[i+4*N];
        nu_con[i] = (theta[i+4*N]+1)/(theta[i+4*N]-2);
    }    
     
    for (i=0; i<N; i++) 
    {
        for (j=1; j<T; j++)
        {
            tmp = (y[j-1]-theta[i])*(y[j-1]-theta[i]);
            tmp = tmp/(h[i]*(theta[i+4*N]-2));
            tmp = 1 + tmp;
            tmp = nu_con[i]/tmp;
            h[i] = theta[3*N+i]*h[i] + theta[i+N] + A[i]*(tmp*(y[j-1]-theta[i])*(y[j-1]-theta[i]) - h[i]);
        }                  
     }
    
    /* Free allocated memory */
    mxFree(A); 
    mxFree(nu_con); 
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mwSignedIndex N, T;                               /* size of matrix */
    double *y, *theta;                                /* input*/
    double *h;                                        /* output */
    
    /* Getting the inputs */
    theta = mxGetPr(prhs[0]);
    y = mxGetPr(prhs[1]);
    
    N = mxGetM(prhs[0]); /* no of parameter draws */
    T = mxGetM(prhs[1]); /* no. of observations */
    
    /* create the output matrices */
    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL); 

    /* get a pointer to the real data in the output matrix */
    h = mxGetPr(plhs[0]);
    
    /* call the function */
    volatility_t_gas_mex(y, N, T, theta, h);
  
}
