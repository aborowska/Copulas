#include "mex.h"
#include "math.h"
#include "matrix.h"

/*******************************************************************  */

void my_norminv(double *ninv, double *y,  mwSignedIndex T, mwSignedIndex D)
{
    mxArray *in_array_ptr, *out_array_ptr; // mxArray * - a pointer to a struct (A POINTER TO A POINTER??)

    in_array_ptr  = mxCreateDoubleMatrix(T, D, mxREAL);      
    
    memcpy(mxGetPr(in_array_ptr), y, T*D*sizeof(double)); // start copying at the double returned by mxGetPr(array_ptr)
    
    mexCallMATLAB(1, &out_array_ptr, 1, &in_array_ptr, "norminv"); // & turns a value into a pointer --> call a MATLAB function with pointer to real matrices
    memcpy(ninv, mxGetPr(out_array_ptr), T*D*sizeof(double)); // start copying at the double returned by d=mxGetPr(plhs[0])
 
    mxDestroyArray(in_array_ptr);
    mxDestroyArray(out_array_ptr);
}

/*******************************************************************  */

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mwSignedIndex N, k, T, D;                           /* size of matrix */
    double *y;                /* input*/
    double *ninv;                                   /* output */

    /* Getting the inputs */
    y = mxGetPr(prhs[0]);
      
    T = mxGetM(prhs[0]); /* no. of observations */
    D = mxGetN(prhs[0]); /* no. of series */
    
    
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(T,D,mxREAL); /* log posterior */
  
    /* get a pointer to the real data in the output matrix */
    ninv = mxGetPr(plhs[0]);
     
    /* call the function */
    my_norminv(ninv, y, T, D);
}
