#include "mex.h"
#include "math.h"
#include "matrix.h"

#define  PI  3.1415926535897932;

void duvt_garch(double x, double mu, double sigma, double df, double *GamMat, mwSignedIndex G, double *pdf)
{
    double c0, c1, c2, c, e, tmp, etmp, df5;
    int ind;
    
    c0 = df+1;
    
    if (c0 <= 100)
    {
        ind = floor(c0*50000) - 1;
        c1 = GamMat[ind];
    }
    else
    {
        c1 = GamMat[G-1];
    }
    
    if (df <= 100)
    {
        df5 = df*50000;
        if ((df5 - floor(df5)) < (floor(df5+1) - df5))
        {
            ind = floor(df5);
        }
        else
        {
            ind = floor(df5 + 1);
        }
        ind = ind -  1;
        c2 = GamMat[ind];
    }
    else
    {
        c2 = GamMat[G-1];
    } 

    c = df*PI;
    c = pow(c,0.5);
    c2 = c2*c;
    c2 = c2*pow(sigma,0.5);
    c = c1/c2;
    e = -0.5*c0; 
    
    tmp = (x-mu)*(x-mu)/sigma;
    tmp = 1 + tmp/df;
    etmp = pow(tmp,e); 
    pdf[0] = exp(log(c) + log(etmp));
}

void loglik_t_garch_noS_hyper_init_mex(double *y, mwSignedIndex N, mwSignedIndex T, double *S,
        double *theta, double *hyper, double *GamMat, mwSignedIndex G, double *d,
        double *Td)
{
    double *rho;
    double h, *pdf; 
    double rhoh;
    mwSignedIndex i, j;
        
    /* Variable size arrays */ 
    rho = mxMalloc((N)*sizeof(double));
    pdf = mxMalloc((1)*sizeof(double));
    
    Td[0] = (double)(T);
    
    /* Initialise */
    for (i=0; i<N; i++)
    {
        rho[i] = (theta[i+4*N]-2)/theta[i+4*N];
    }
         
    /* loglik */
    for (i=0; i<N; i++) 
    {         
        d[i] = 0;
        h = S[0]; 
        for (j=1; j<T; j++)
        {   
            h = theta[2*N+i]*h + theta[i] + theta[N+i]*(y[j-1]-theta[3*N+i])*(y[j-1]-theta[3*N+i]);   
            rhoh = rho[i]*h;
            duvt_garch(y[j], theta[i+3*N], rhoh, theta[i+4*N], GamMat, G, pdf);
            d[i] = d[i] + log(pdf[0]);
        }      
        d[i] = -d[i]/T;
     }
    
    /* Free allocated memory */
    mxFree(rho); 
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mwSignedIndex N, T, G;                              /* size of matrix */
    double *y, *S, *theta, *hyper;                      /* input*/
    double *GamMat;                                     /* input */
    double *d;                                          /* output */
    double *Td;

    /* Getting the inputs */
    theta = mxGetPr(prhs[0]);
    y = mxGetPr(prhs[1]);
    S = mxGetPr(prhs[2]);
    GamMat = mxGetPr(prhs[3]);
    hyper = mxGetPr(prhs[4]); /* hyperparameter on degrees of freedom*/
    
    N = mxGetM(prhs[0]); /* no of parameter draws */
    T = mxGetM(prhs[1]); /* no. of observations */
    G = mxGetM(prhs[3]);
        
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL); 
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);

    /* get a pointer to the real data in the output matrix */
    d = mxGetPr(plhs[0]);
    Td = mxGetPr(plhs[1]);
    
    /* call the function */
    loglik_t_garch_noS_hyper_init_mex(y, N, T, S, theta, hyper, GamMat, G, d, Td); 
}