#include "mex.h"
#include "math.h"
#include "matrix.h"

#define PI      3.14159265358979323846
#define Log2PI  1.83787706640934548356
#define LogPI   1.14472988584940017414

/*******************************************************************  */

void my_tinv(double *tinv, double *y, double *nu, mwSignedIndex T, mwSignedIndex D)
{
    mxArray *in_array_ptr[2], *out_array_ptr; // mxArray * - a pointer to a struct (A POINTER TO A POINTER??)

    in_array_ptr[0] = mxCreateDoubleMatrix(T, D, mxREAL);      
    in_array_ptr[1] = mxCreateDoubleMatrix(1, 1, mxREAL);  
    
    memcpy(mxGetPr(in_array_ptr[0]), y, T*D*sizeof(double)); // start copying at the double returned by mxGetPr(array_ptr)
    memcpy(mxGetPr(in_array_ptr[1]), nu, sizeof(double)); // start copying at the double returned by mxGetPr(array_ptr)
    
    mexCallMATLAB(1, &out_array_ptr, 2, in_array_ptr, "tinv"); // & turns a value into a pointer --> call a MATLAB function with pointer to real matrices
    memcpy(tinv, mxGetPr(out_array_ptr), T*D*sizeof(double)); // start copying at the double returned by d=mxGetPr(plhs[0])
 
    mxDestroyArray(in_array_ptr);
    mxDestroyArray(out_array_ptr);
}

/*******************************************************************  */

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mwSignedIndex N, k, T, D;                           /* size of matrix */
    double *y, *nu;                /* input*/
    double *tinv;                                   /* output */

    /* Getting the inputs */
    nu = mxGetPr(prhs[0]);
    y = mxGetPr(prhs[1]);
      
    N = mxGetM(prhs[0]); /* no of parameter draws */
    T = mxGetM(prhs[1]); /* no. of observations */
    D = mxGetN(prhs[1]); /* no. of series */
    
    
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(T,D,mxREAL); /* log posterior */
  
    /* get a pointer to the real data in the output matrix */
    tinv = mxGetPr(plhs[0]);
     
    /* call the function */
    my_tinv(tinv, y, nu, T, D);
}
