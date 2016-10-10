#include "mex.h"
#include "math.h"
#include "matrix.h"

#define  PI  3.1415926535897932;
// const double M = -1e100;

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
     /*   c1 = tgamma(0.5*c0); */
        c1 = GamMat[G-1];
    }
    
    if (df <= 100)
    {
//         ind = floor(df[0]*50000) - 1; /* useround insteadas follows: */
        df5 = df*50000;
        if ((df5 - floor(df5)) < (floor(df5+1) - df5))
        {
            ind = floor(df5);
        }
        else
        {
            ind = floor(df5 + 1);
        }
//     	mexPrintf("ind df = %i\n",ind);       
        ind = ind -  1;
        c2 = GamMat[ind];
    }
    else
    {
        c2 = GamMat[G-1];
    } 

    c = df*PI;
    c = pow(c,0.5);
//     mexPrintf("c = %6.4f\n",c);    
    c2 = c2*c;
    c2 = c2*pow(sigma,0.5);
    c = c1/c2;
    e = -0.5*c0; 
    
    tmp = (x-mu)*(x-mu)/sigma;
    tmp = 1 + tmp/df;
    etmp = pow(tmp,e); 
    pdf[0] = exp(log(c) + log(etmp));
}

void prior_t_garch_hyper(double *theta, double *hyper, 
        mwSignedIndex N, mwSignedIndex *r1, double *r2)
{
    mwSignedIndex i;
    
    /* Variable size arrays */
    /*  c1 = malloc((N)*sizeof(double));             */
    
    for (i=0; i<N; i++)
    {
        r1[i] = 1;
        if (theta[i] <= 0) // omega>0
        {
            r1[i] = 0;
        }   
        if ((theta[i+N] <0 ) || (theta[i+N] >= 1)) // 0<=alpha<1 
        {
            r1[i] = 0;
        }
        if ((theta[i+2*N] < 0) || (theta[i+2*N] >= 1)) // 0<=beta<1
        {
            r1[i] = 0;
        }
        if (theta[i+N] + theta[i+2*N] >= 1) //alpha+beta<1
        {
            r1[i] = 0;
        }
        if (theta[i+4*N] <= 2) //nu>2
        {
            r1[i] = 0;
        }     
        if (r1[i] == 1)
        {
            r2[i] = log(hyper[0]) - hyper[0]*(theta[i+4*N] - 2);
        }
    }
}

void posterior_t_garch_noS_hyper_init_mex(double *y, mwSignedIndex N, mwSignedIndex T, double *S,
        double *theta, double *hyper, double *GamMat, mwSignedIndex G, double *d,
        double *Td)
{
    mwSignedIndex *r1;
    double *r2; 
//    double *omega, *rho;
    double *rho;
    double h, *pdf; 
    double rhoh;
    mwSignedIndex i, j;
        
    /* Variable size arrays */
    r1 = mxMalloc((N)*sizeof(mwSignedIndex));              
    r2 = mxMalloc((N)*sizeof(double));              
 //   omega = mxMalloc((N)*sizeof(double));   
    rho = mxMalloc((N)*sizeof(double));
    pdf = mxMalloc((1)*sizeof(double));
    
    Td[0] = (double)(T);
    
    prior_t_garch_hyper(theta, hyper, N, r1, r2);

    /* Initialise */
    for (i=0; i<N; i++)
    {
//        omega[i] = S[0]*(1-theta[i]-theta[N+i]);  
        rho[i] = (theta[i+4*N]-2)/theta[i+4*N];
//         mexPrintf("omega[%i] = %6.4f\n",i,omega[i] );       
//         mexPrintf("rho[%i] = %6.4f\n",i,rho[i]);  
//         mexPrintf("alpha[%i] = %6.4f\n", i, theta[i]); 
//         mexPrintf("beta[%i] = %6.4f\n", i, theta[i+N]);
//         mexPrintf("mu[%i] = %6.4f\n", i, theta[2*N+i]);
//         mexPrintf("nu[%i] = %6.4f\n", i, theta[i+3*N]);
    }

         
    /* PDF */
    for (i=0; i<N; i++) 
    {         
        
        if (r1[i]==1)
        {
            d[i] = r2[i];
            h = S[0]; 
//             mexPrintf("d[%i] = %6.4f\n",i,d[i]);       
//             mexPrintf("h = %6.4f\n", h);   
            for (j=1; j<T; j++)
            {   
//                 mexPrintf("i = %i\n",i); 
//                 mexPrintf("j = %i\n",j); 
//                 mexPrintf("y[%i] = %6.4f\n", j-1, y[j-1]);  
//                 
//                 mexPrintf("alpha[%i] = %6.4f\n", i, theta[i]); 
//                 mexPrintf("beta[%i] = %6.4f\n", i, theta[i+N]);
//                 mexPrintf("mu[%i] = %6.4f\n", i, theta[2*N+i]);
//                 mexPrintf("nu[%i] = %6.4f\n", i, theta[i+3*N]);
// //     /*            h[i] = beta[i]*h[i] + omega[i] + alpha[i]*(y[j-1]-mu[i])*(y[j-1]-mu[i]);     */    
//                h = theta[N+i]*h + omega[i] + theta[i]*(y[j-1]-theta[2*N+i])*(y[j-1]-theta[2*N+i]);   
                h = theta[2*N+i]*h + theta[i] + theta[N+i]*(y[j-1]-theta[3*N+i])*(y[j-1]-theta[3*N+i]);   
//                 mexPrintf("h = %6.4f\n", h);   
                rhoh = rho[i]*h;
//                 mexPrintf("rhoh = %6.4f\n", rhoh);  
                duvt_garch(y[j], theta[i+3*N], rhoh, theta[i+4*N], GamMat, G, pdf);
//                 mexPrintf("pdf[%i] = %16.14f\n", j, pdf[0]);  
                d[i] = d[i] + log(pdf[0]);
            }      
            d[i] = -d[i]/T;
        }
        else
        {
//             d[i] = M;  
            d[i] = mxGetInf();
        }    
//         mexPrintf("d[%i] = %6.4f\n",i,d[i]);  
     }
    
    /* Free allocated memory */
    mxFree(r1); 
    mxFree(r2); 
    mxFree(rho); 
//    mxFree(omega);
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
    
//     mexPrintf("N = %i\n", N);  
//     mexPrintf("T = %i\n", T); 
//     mexPrintf("G = %i\n", G); 
        
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL); 
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);

    /* get a pointer to the real data in the output matrix */
    d = mxGetPr(plhs[0]);
    Td = mxGetPr(plhs[1]);
    
    /* call the function */
    posterior_t_garch_noS_hyper_init_mex(y, N, T, S, theta, hyper, GamMat, G, d, Td);
  
}
