#include "mex.h"
#include "math.h"
#include "matrix.h"

#define  PI     3.1415926535897932;
#define Log2PI  1.837877066409345;
/******************************************************************* */

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

/*******************************************************************  */


void prior_gen_gas(double *theta, double *hyper, 
        mwSignedIndex N, mwSignedIndex k, mwSignedIndex *r1, double *r2)
{
    mwSignedIndex i;
    
    /* Variable size arrays */
    
    for (i=0; i<N; i++)
    {
        r1[i] = 1;
        r2[i] = -mxGetInf();
       
        if (theta[i+N] <= 0) // omega>0
        {
            r1[i] = 0;
        }   
        if ((theta[i+3*N] < 0) || (theta[i+3*N] >= 1)) // 0<=B<1
        {
            r1[i] = 0;
        }
        if (k == 5)
        {
            if (theta[i+4*N] <= 2) //nu>2
            {
                r1[i] = 0;
            }     
            if (r1[i] == 1)
            {
                r2[i] = log(hyper[0]) - hyper[0]*(theta[i+4*N] - 2);
            }        
        }
        else
        {
            if (r1[i] == 1)
            {
                r2[i] = 1;
            }        
        }

    }
}

/*******************************************************************  */


void posterior_gen_gas_mex(double *y, mwSignedIndex N, mwSignedIndex k, mwSignedIndex T,
        double *theta, double *hyper, double *link, double *scale, double *GamMat, mwSignedIndex G, 
        double *d, double *f)
{
    mwSignedIndex *r1;
    double *r2; 
    double *rho, *nu_con;
    double sigma2, *pdf; 
    double rhos, tmp, S, scaled_score;
    mwSignedIndex i, j;
        
    /* Variable size arrays */
    r1 = mxMalloc((N)*sizeof(mwSignedIndex));              
    r2 = mxMalloc((N)*sizeof(double));            
    rho = mxMalloc((N)*sizeof(double));
    nu_con = mxMalloc((N)*sizeof(double));

    pdf = mxMalloc((1)*sizeof(double));
    
    prior_gen_gas(theta, hyper, N, k, r1, r2);

    /* Initialise */
    if (k == 5)
    {
        for (i=0; i<N; i++)
        {
            rho[i] = (theta[i+4*N]-2)/theta[i+4*N];
            nu_con[i] = 2*(theta[i+4*N]+3)/theta[i+4*N];
        }
    }

      
    /* PDF */
    for (i=0; i<N; i++) 
    {  
        if (r1[i]==1)
        {
            d[i] = r2[i];
            f[i] = theta[i+N]/(1-theta[i+3*N]);
            if (link[0]) /* linear link*/
            {
                sigma2 = f[i];             
            }
            else    /* exp link*/
            {
                sigma2 = exp(f[i]);             
            }

            if (k == 5)
            {
                rhos = rho[i]*sigma2;
                duvt_garch(y[0], theta[i], rhos, theta[i+4*N], GamMat, G, pdf);   
                pdf[0] = log(pdf[0]);
            }
            else
            {
                pdf[0] =  (y[0]-theta[i])*(y[0]-theta[i])/sigma2;   
                pdf[0] = Log2PI + log(sigma2) + pdf[0];   
                pdf[0] = -0.5*pdf[0];   
            }
            d[i] = d[i] + pdf[0];
            
            
            for (j=1; j<T; j++)
            {   
                if (k == 5) 
                {
                    S = nu_con[i]*sigma2*sigma2; /* S == inv_fisher */                    
                            
                    tmp = (y[j-1]-theta[i])*(y[j-1]-theta[i]);
                    tmp = tmp/sigma2;
                    tmp = tmp + theta[i+4*N]-2;             
                    tmp = (theta[i+4*N]+1)/tmp;   
                    
                    scaled_score = (tmp*(y[j-1] - theta[i])*(y[j-1] - theta[i])/sigma2 - 1);
                    scaled_score = scaled_score/(2*sigma2);                           
                }               
                else
                {
                    S = 2*sigma2*sigma2; /* S == inv_fisher */
                    scaled_score = ((y[j-1] - theta[i])*(y[j-1] - theta[i])/sigma2 - 1);
                    scaled_score = scaled_score/(2*sigma2);                            
                }

                if (!link[0]) /* apply chain rule for the exp link */
                {
//                     mexPrintf("chain rule\n");                     
                    S = S/(sigma2*sigma2);    
                    scaled_score = scaled_score*sigma2;                        
                }
                
                if (!scale[0]) /* apply sqrt scaling*/
                {
//                    mexPrintf("sqrt scaling\n"); 
                   S = pow(S,0.5);
                }  
                
                scaled_score = S*scaled_score;
                f[N*j + i] = theta[i+N] + theta[2*N+i]*scaled_score + theta[3*N+i]*f[N*(j-1) + i]  ; 
                if (link[0]) /* linear link*/
                {
//                     mexPrintf("linear link\n");                     
                    sigma2 = f[N*j + i];             
                }
                else    /* exp link*/
                {
//                     mexPrintf("exp link\n"); 
                    sigma2 = exp(f[N*j + i]);             
                } 
                
                if (k == 5)
                {
//                     mexPrintf("student\n");                     
                    rhos = rho[i]*sigma2;
                    duvt_garch(y[j] , theta[i], rhos, theta[i+4*N], GamMat, G, pdf);   
                    pdf[0] = log(pdf[0]);
                }
                else
                {
//                     mexPrintf("normal\n");                     
                    pdf[0] = (y[j]-theta[i])*(y[j]-theta[i])/sigma2;                       
                    pdf[0] = Log2PI + log(sigma2) + pdf[0];
                    pdf[0] = -0.5*pdf[0];   
                    
                }     
                d[i] = d[i] + pdf[0];
            }         
        }
        else
        {
            d[i] = -mxGetInf();
        }    
     }
    
    /* Free allocated memory */
    mxFree(r1); 
    mxFree(r2); 
    mxFree(rho); 
    mxFree(nu_con); 
    mxFree(pdf);
}

/*******************************************************************  */

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mwSignedIndex N, k, T, G;                           /* size of matrix */
    double *y, *theta, *hyper;                          /* input*/
    double *link, *scale;                               /* input*/
    double *GamMat;                                     /* input */
    double *d, *f;                                      /* output */

    /* Getting the inputs */
    theta = mxGetPr(prhs[0]);
    y = mxGetPr(prhs[1]);
    hyper = mxGetPr(prhs[2]); /* hyperparameter on degrees of freedom if student's t errors */
    link = mxGetPr(prhs[3]);
    scale = mxGetPr(prhs[4]);
    GamMat = mxGetPr(prhs[5]);
    
    N = mxGetM(prhs[0]); /* no of parameter draws */
    k = mxGetN(prhs[0]); /* draw dimension */    
    T = mxGetM(prhs[1]); /* no. of observations */
    G = mxGetM(prhs[5]);
    
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL); /* log posterior */
    plhs[1] = mxCreateDoubleMatrix(N,T,mxREAL); /* corresponding signals */

    /* get a pointer to the real data in the output matrix */
    d = mxGetPr(plhs[0]);
    f = mxGetPr(plhs[1]);
    
    /* call the function */
    posterior_gen_gas_mex(y, N, k, T, theta, hyper, link, scale, GamMat, G, d, f);
}
