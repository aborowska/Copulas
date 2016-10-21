#include "mex.h"
#include "math.h"
#include "matrix.h"

#define  PI     3.14159265358979323846
#define Log2PI  1.83787706640934548356
#define LogPI   1.14472988584940017414
/*******************************************************************  */

void my_gamma(double *gam, mwSignedIndex N)
{
    mxArray *in_array_ptr, *out_array_ptr; // mxArray * - a pointer to a struct (A POINTER TO A POINTER??)
   
    in_array_ptr = mxCreateDoubleMatrix(2*N, 1, mxREAL);  
    
    memcpy(mxGetPr(in_array_ptr), gam, 2*N*sizeof(double)); // start copying at the double returned by mxGetPr(array_ptr)
    mexCallMATLAB(1, &out_array_ptr, 1, &in_array_ptr, "gamma"); // & turns a value into a pointer --> call a MATLAB function with pointer to real matrices
    memcpy(gam, mxGetPr(out_array_ptr), 2*N*sizeof(double)); // start copying at the double returned by d=mxGetPr(plhs[0])
 
    mxDestroyArray(in_array_ptr);
    mxDestroyArray(out_array_ptr);
}

/******************************************************************* */
void loglik_gen_gas_mex(double *y, mwSignedIndex N, mwSignedIndex k, mwSignedIndex T,
        double *theta, double *link, double *scale,
        double *d, double *f)
{
 
    double *nu_con, *gam, *gam2;
    double y2, sigma2, pdf; 
    double tmp, S, scaled_score;
    mwSignedIndex i, j;
    
    /* Variable size arrays */
    if (k == 5)
    {
        nu_con = mxMalloc((N)*sizeof(double));
        gam = mxMalloc(2*(N)*sizeof(double));
    }
    
    /* Initialise */
    if (k == 5)
    {
        for (i=0; i<N; i++)
        {
            gam[2*i] = (theta[i+4*N]+1)/2;
            gam[2*i+1] = theta[i+4*N]/2;            
            nu_con[i] = 2*(theta[i+4*N]+3)/theta[i+4*N];
        }
        my_gamma(gam,N);
    }
    
    /* Get the constants for Student's distribution */
     
    /* PDF */
    for (i=0; i<N; i++) 
    {  
        pdf = 0;
        f[i] = theta[i+N]/(1-theta[i+3*N]);
        if (link[0]) /* linear link*/
        {
            sigma2 = f[i];             
        }
        else    /* exp link*/
        {
            sigma2 = exp(f[i]);             
        }
        
        y2 = (y[0]-theta[i])*(y[0]-theta[i])/sigma2;
        if (k == 5)
        {
            y2 = y2/(theta[i+4*N] - 2);
            pdf = pdf + log(gam[2*i]);           
            pdf = pdf - log(gam[2*i+1]);          
        }
        else
        {
            pdf = -0.5*(Log2PI + log(sigma2) + y2);   
        }
        d[i] = pdf;


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
            
            y2 = (y[j]-theta[i])*(y[j]-theta[i])/sigma2;
            if (k == 5)
            {           
//                     mexPrintf("student\n");
                y2 = y2/(theta[i+4*N] - 2);
                pdf = - 0.5*(log(theta[i+4*N]-2) + log(PI) + log(sigma2) + (theta[i+4*N]+1)*log(1+y2));
                pdf = pdf + log(gam[2*i]);           
                pdf = pdf - log(gam[2*i+1]); 
            }
            else
            {
//                     mexPrintf("normal\n");                     
                pdf = -0.5*(Log2PI + log(sigma2) + y2);   
            }     
            d[i] = d[i] + pdf;
        }     
        d[i] = -d[i]/T;    
     }
    
    /* Free allocated memory */
    if (k == 5)
    {    
        mxFree(nu_con); 
        mxFree(gam); 
    }
}

/*******************************************************************  */

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mwSignedIndex N, k, T;                           /* size of matrix */
    double *y, *theta, *link, *scale;                /* input*/
    double *d, *f;                                   /* output */

    /* Getting the inputs */
    theta = mxGetPr(prhs[0]);
    y = mxGetPr(prhs[1]);
    link = mxGetPr(prhs[2]);
    scale = mxGetPr(prhs[3]);
      
    N = mxGetM(prhs[0]); /* no of parameter draws */
    k = mxGetN(prhs[0]); /* draw dimension */    
    T = mxGetM(prhs[1]); /* no. of observations */
      
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL); /* log posterior */
    plhs[1] = mxCreateDoubleMatrix(N,T,mxREAL); /* corresponding signals */

    /* get a pointer to the real data in the output matrix */
    d = mxGetPr(plhs[0]);
    f = mxGetPr(plhs[1]);
    
    /* call the function */
    loglik_gen_gas_mex(y, N, k, T, theta, link, scale, d, f);
}
