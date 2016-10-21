#include "mex.h"
#include "math.h"
#include "matrix.h"

/*******************************************************************  */

void my_norminv(double *ninv, double *u, mwSignedIndex T)
{
    mxArray *in_array_ptr, *out_array_ptr; // mxArray * - a pointer to a struct (A POINTER TO A POINTER??)

    in_array_ptr = mxCreateDoubleMatrix(T, 2, mxREAL);      
    
    memcpy(mxGetPr(in_array_ptr), u, T*2*sizeof(double)); // start copying at the double returned by mxGetPr(array_ptr)
    
    mexCallMATLAB(1, &out_array_ptr, 1, &in_array_ptr, "norminv"); // & turns a value into a pointer --> call a MATLAB function with pointer to real matrices
    memcpy(ninv, mxGetPr(out_array_ptr), T*2*sizeof(double)); // start copying at the double returned by d=mxGetPr(plhs[0])
 
    mxDestroyArray(in_array_ptr);
    mxDestroyArray(out_array_ptr);
}

/*******************************************************************  */

void my_tinv(double *tinv, double *u, double *nu, mwSignedIndex T)
{
    mxArray *in_array_ptr[2], *out_array_ptr; // mxArray * - a pointer to a struct (A POINTER TO A POINTER??)

    in_array_ptr[0] = mxCreateDoubleMatrix(T, 2, mxREAL);      
    in_array_ptr[1] = mxCreateDoubleMatrix(1, 1, mxREAL);  
    
    memcpy(mxGetPr(in_array_ptr[0]), u, T*2*sizeof(double)); // start copying at the double returned by mxGetPr(array_ptr)
    memcpy(mxGetPr(in_array_ptr[1]), nu, sizeof(double)); // start copying at the double returned by mxGetPr(array_ptr)
    
    mexCallMATLAB(1, &out_array_ptr, 2, in_array_ptr, "tinv"); // & turns a value into a pointer --> call a MATLAB function with pointer to real matrices
    memcpy(tinv, mxGetPr(out_array_ptr), T*2*sizeof(double)); // start copying at the double returned by d=mxGetPr(plhs[0])
 
    mxDestroyArray(in_array_ptr);
    mxDestroyArray(out_array_ptr);
}

/*******************************************************************  */

void my_gamma(double *gam, mwSignedIndex N)
{
    mxArray *in_array_ptr, *out_array_ptr; // mxArray * - a pointer to a struct (A POINTER TO A POINTER??)
   
    in_array_ptr = mxCreateDoubleMatrix(3*N, 1, mxREAL);  
    
    memcpy(mxGetPr(in_array_ptr), gam, 3*N*sizeof(double)); // start copying at the double returned by mxGetPr(array_ptr)
    mexCallMATLAB(1, &out_array_ptr, 1, &in_array_ptr, "gamma"); // & turns a value into a pointer --> call a MATLAB function with pointer to real matrices
    memcpy(gam, mxGetPr(out_array_ptr), 3*N*sizeof(double)); // start copying at the double returned by d=mxGetPr(plhs[0])
 
    mxDestroyArray(in_array_ptr);
    mxDestroyArray(out_array_ptr);
}

/******************************************************************* */
void loglik_gen_copula_gas_mex(double *y, mwSignedIndex N, mwSignedIndex k, mwSignedIndex T,
        double *theta, double *link, double *scale,
        double *d, double *f, double *rho)
{
 
    double *z, *gam;
    double pdf; /*,rho;*/
    double tmp, rr, S, scaled_score;
    mwSignedIndex i, j;
    
    /* Variable size arrays */
    z = mxMalloc(2*(T)*sizeof(double));    
    if (k == 4)
    {        
        gam = mxMalloc(3*(N)*sizeof(double));
    }
    
    /* Initialise */
    if (k == 4)   /* Get the constants for Student's distribution */
    {
        for (i=0; i<N; i++)
        {        
            gam[3*i] = (theta[i+3*N]+2)/2;
            gam[3*i+1] = theta[i+3*N]/2;            
            gam[3*i+2] = (theta[i+3*N]+1)/2;   
        }
        my_gamma(gam,N);
    }
    else
    {
        my_norminv(z, y, T);
    }
    
     
    /* PDF */
    for (i=0; i<N; i++) 
    {  
        if (k == 4)
        {
            my_tinv(z, y, &theta[i+3*N], T);          
        }
        
        pdf = 0;
        f[i] = theta[i]/(1-theta[i+2*N]);
        
        if (link[0]) /* KLS link*/
        {
            rho[i] = (1 - exp(-f[i]))/(1 + exp(-f[i]));    
//             mexPrintf("KLS link\n");
        }
        else    /* NAIS link*/
        {
            rho[i] = 1/(1 + exp(-f[i]));     
//             mexPrintf("NAIS link\n");            
        }
        
        if (k == 4)
        {
            pdf = z[0]*z[0] + z[0+T]*z[0+T] - 2*rho[i]*z[0]*z[0+T];
            pdf = pdf/(theta[3*N+i]*(1-rho[i]*rho[i]));
            pdf = (theta[3*N+i] + 2)*log(1 + pdf);
            pdf = pdf + log(1 - rho[i]*rho[i]) - (theta[3*N+i] + 1)*(log(1 + z[0]*z[0]/theta[3*N+i]) + log(1 + z[0+T]*z[0+T]/theta[3*N+i]));  
            pdf = -0.5*pdf;                   
            pdf = pdf + log(gam[3*i]);           
            pdf = pdf + log(gam[3*i+1]); 
            pdf = pdf - 2*log(gam[3*i+2]);                
//             mexPrintf("z[%i] = %10.8f\n",1,z[0]);  
//             mexPrintf("z[%i] = %10.8f\n",2,z[T]); 
//             mexPrintf("rho[%i] = %10.8f\n",0,rho[i]);  
//             mexPrintf("pdf[%i] = %10.8f\n",0,pdf);  
            
        }
        else
        {
            pdf = z[0]*z[0] + z[0+T]*z[0+T] - 2*rho[i]*z[0]*z[0+T];
            pdf = pdf/(1-rho[i]*rho[i]);
            pdf = pdf + log(1-rho[i]*rho[i]) - z[0]*z[0] - z[0+T]*z[0+T]; 
            pdf = -0.5*pdf;   
        }
        d[i] = pdf;


        for (j=1; j<T; j++)
        {   
            if (k == 4) 
            {
                rr = rho[N*(j-1) + i];
                S = (theta[3*N+i]+4)*(1 - rr*rr)*(1 - rr*rr)/(theta[3*N+i] + 2 + theta[3*N+i]*rr*rr); /* S == inv_fisher */
                tmp = z[j-1]*z[j-1] + z[j-1+T]*z[j-1+T] - 2*rr*z[j-1]*z[j-1+T];
                tmp = theta[3*N+i] + tmp/(1-rr*rr);
                tmp = (theta[i+3*N]+2)/tmp;   
             
                scaled_score = (1 + rr*rr)*(tmp*z[j-1]*z[j-1+T] - rr) - rr*(tmp*z[j-1]*z[j-1] + tmp*z[j-1+T]*z[j-1+T] - 2);    
                scaled_score = scaled_score/((1 - rr*rr)*(1 - rr*rr));                            
            }               
            else
            {
                S = (1-rho[N*(j-1) + i]*rho[N*(j-1) + i])*(1-rho[N*(j-1) + i]*rho[N*(j-1) + i])/(1+rho[N*(j-1) + i]*rho[N*(j-1) + i]); /* S == inv_fisher */
                scaled_score = (1+rho[N*(j-1) + i]*rho[N*(j-1) + i])*(z[j-1]*z[j-1+T]-rho[N*(j-1) + i]) - rho[N*(j-1) + i]*(z[j-1]*z[j-1] + z[j-1+T]*z[j-1+T] - 2);    
                scaled_score = scaled_score/((1-rho[N*(j-1) + i]*rho[N*(j-1) + i])*(1-rho[N*(j-1) + i]*rho[N*(j-1) + i]));                
            }

            if (link[0]) /* apply chain rule */
            {
               // mexPrintf("KLS   link\n");                            
                tmp = 2/(2 + exp(f[N*(j-1) + i]) + exp(-f[N*(j-1) + i]));                        
            }
            else
            {
//                 mexPrintf("NAIS link\n"); 
                tmp = 1/(2 + exp(f[N*(j-1) + i]) + exp(-f[N*(j-1) + i]));                 
            }
            S = S/(tmp*tmp);    
            scaled_score = tmp*scaled_score;
 
            if (!scale[0]) /* apply sqrt scaling*/
            {
//                 mexPrintf("sqrt scaling\n");            
                S = pow(S,0.5);
            }  

            scaled_score = S*scaled_score;
            f[N*j + i] = theta[i] + theta[N+i]*scaled_score + theta[2*N+i]*f[N*(j-1) + i]; 
            
            if (link[0]) /* KLS link*/
            {
//                 mexPrintf("KLS link\n");                        
                rho[N*j+ i] = (1 - exp(-f[N*j + i]))/(1 + exp(-f[N*j + i]));             
            }
            else    /* NAIS link*/
            {
//                 mexPrintf("NAIS link\n");                            
                rho[N*j+ i] = 1/(1 + exp(-f[N*j + i]));                  
            }
            
            if (k == 4)
            {           
//                     mexPrintf("student\n");
                pdf = z[j]*z[j] + z[j+T]*z[j+T] - 2*rho[N*j+ i]*z[j]*z[j+T];
                pdf = pdf/(theta[3*N+i]*(1-rho[N*j+ i]*rho[N*j+ i]));
                pdf = (theta[3*N+i] + 2)*log(1 + pdf);
                pdf = pdf + log(1 - rho[N*j+i]*rho[N*j+i]) - (theta[3*N+i] + 1)*(log(1 + z[j]*z[j]/theta[3*N+i]) + log(1 + z[j+T]*z[j+T]/theta[3*N+i]));  
                pdf = -0.5*pdf;                   
                pdf = pdf + log(gam[3*i]);           
                pdf = pdf + log(gam[3*i+1]); 
                pdf = pdf - 2*log(gam[3*i+2]);         
//                 mexPrintf("pdf[%i] = %10.8f\n",j,pdf);                                           
            }
            else
            {
//                     mexPrintf("normal\n");                     
                pdf = z[j]*z[j] + z[j+T]*z[j+T] - 2*rho[N*j+ i]*z[j]*z[j+T];
                pdf = pdf/(1-rho[N*j+ i]*rho[N*j+ i]);
                pdf = pdf + log(1-rho[N*j+ i]*rho[N*j+ i]) - z[j]*z[j] - z[j+T]*z[j+T]; 
                pdf = -0.5*pdf;     
            }     
            d[i] = d[i] + pdf;
        }     
        d[i] = -d[i]/T;    
     }
    
    /* Free allocated memory */
    mxFree(z);     
    if (k == 4)
    {    
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
    double *d, *f, *rho;                                   /* output */

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
    plhs[2] = mxCreateDoubleMatrix(N,T,mxREAL); /* corresponding signals */

    /* get a pointer to the real data in the output matrix */
    d = mxGetPr(plhs[0]);
    f = mxGetPr(plhs[1]);
    rho = mxGetPr(plhs[2]);
                
    /* call the function */
    loglik_gen_copula_gas_mex(y, N, k, T, theta, link, scale, d, f, rho);
}
