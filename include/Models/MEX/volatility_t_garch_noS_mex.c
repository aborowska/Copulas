#include "mex.h"
#include "math.h"

void volatility_t_garch(double *y, mwSignedIndex N, mwSignedIndex T, double *S,
        double *theta, double *h)
{
//    double *omega; 
    mwSignedIndex i, j;
        
    /* Variable size arrays */
//    omega = malloc((N)*sizeof(double));              
    
    /* Initialise */
    for (i=0; i<N; i++)
    {
/*        omega[i] = S[0]*(1-alpha[i]-beta[i]); */
//        omega[i] = S[0]*(1-theta[i]-theta[N+i]);
        h[i] = S[0]; 
    }    
     
    for (i=0; i<N; i++) 
    {
        for (j=1; j<T; j++)
        {   
/*            h[i] = beta[i]*h[i] + omega[i] + alpha[i]*(y[j-1]-mu[i])*(y[j-1]-mu[i]);     */    
//            h[i] = theta[N+i]*h[i] + omega[i] + theta[i]*(y[j-1]-theta[2*N+i])*(y[j-1]-theta[2*N+i]);         
            h[i] = theta[2*N+i]*h[i] + theta[i] + theta[N+i]*(y[j-1]-theta[3*N+i])*(y[j-1]-theta[3*N+i]);         
        }                  
     }
    
    /* Free allocated memory */
//    free(omega);  
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mwSignedIndex N, T;                               /* size of matrix */
    double *y, *S, *theta;                            /* input*/
    double *h;                                        /* output */
    
    /* Getting the inputs */
    theta = mxGetPr(prhs[0]);
    y = mxGetPr(prhs[1]);
    S = mxGetPr(prhs[2]);
    
    N = mxGetM(prhs[0]); /* no of parameter draws */
    T = mxGetM(prhs[1]); /* no. of observations */
    
    /* create the output matrices */
    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL); 

    /* get a pointer to the real data in the output matrix */
    h = mxGetPr(plhs[0]);
    
    /* call the function */
    volatility_t_garch(y, N, T, S, theta, h);
  
}
