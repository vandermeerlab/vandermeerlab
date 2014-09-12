

#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>



void Usage(void) ;


/******************************************************************************
  INTERFACE FUNCTION
  */

void mexFunction(
   int           nlhs,           /* number of expected outputs */
   mxArray       *plhs[],        /* array of pointers to output arguments */
   int           nrhs,           /* number of inputs */
   const mxArray *prhs[])         /* array of pointers to input arguments */
{

      double 			  *dblptr;
      int               i,j,k,z;
      short int		  tempmax;
      double      *values;
      double      *lookupvalues;        
      double      *directionpr;
      int      direction;
      int			  tmpindex, bound1, bound2;
      double			  *output;
      double                     timeratio;
      int				  nvalues, nlookupvalues, lastnonzero, nextnonzero;

      /* Check numbers of arguments */
      if (!((nrhs == 2)||(nrhs == 3)) || !((nlhs == 1)||(nlhs == 0)) ) {
         Usage();
      }
   

      values = mxGetData(prhs[0]);
      nvalues = mxGetM(prhs[0]) * mxGetN(prhs[0]);

      lookupvalues = mxGetData(prhs[1]);
      nlookupvalues = mxGetM(prhs[1]) * mxGetN(prhs[1]);
      
      if (nrhs == 3) {
         if ((mxGetM(prhs[2]) * mxGetN(prhs[2])) != 1) {
            mexErrMsgTxt("direction is either -1 or 0 or 1 \n");
         }
         directionpr = mxGetData(prhs[2]);
         direction = directionpr[0];
         if ( !((direction == -1)||(direction == 0)||(direction == 1)) ) {
            mexErrMsgTxt("direction is either -1 or 0 or 1\n");
         }
      }
      else {
         direction = 0;
      }
               
      output = (double *) mxCalloc(nvalues, sizeof(double));
      if (nlookupvalues > 0) {     
         for (i = 0; i < nvalues; i++) {
            if ((values[i] < lookupvalues[nlookupvalues-1]) && (values[i] > lookupvalues[0])) {
                  
                  if (nlookupvalues > 2) {
                     tmpindex = (nlookupvalues/2);
                  }
                  else {
                     tmpindex = 0;
                  }
                  
                  bound1 = 0;
                  bound2 = nlookupvalues-1;
                  /*keep dividing the range in half until the two bounds converge around the value*/
                  while ((tmpindex != bound1) && (tmpindex != bound2)) {
                           if (values[i] < lookupvalues[tmpindex]) {
                                 bound2 = tmpindex;
                                 tmpindex = ((tmpindex-bound1)/2) + bound1;
                           }
                           else {
                                 bound1 = tmpindex;
                                 tmpindex = ((bound2-tmpindex)/2)+tmpindex;
                           }
                  }
                  if (direction == 1) {
                     output[i] = tmpindex+2;
                  }
                  else if ((fabs(lookupvalues[tmpindex]-values[i]) < fabs(lookupvalues[tmpindex+1]-values[i])) || (direction == -1) ) {
                     
                     output[i] = tmpindex+1;
                  }
                  else {
                     output[i] = tmpindex+2;
                     
                  }
                                       
            }
            else if (values[i] >= lookupvalues[nlookupvalues-1]){		  
                  output[i] = nlookupvalues;                 
            }
            else if (values[i] <= lookupvalues[0]) {
                  output[i] = 1;   						         
            }
         }   
      }
      else {
         nvalues = 0;
      }
                                                      
      plhs[0] = mxCreateDoubleMatrix(nvalues,1, mxREAL);

      dblptr = mxGetPr(plhs[0]);

      for (i = 0; i < nvalues; i++) {
               *dblptr = output[i];
               dblptr++;
      }
  
      return;
}

void Usage(void)
{
   mexErrMsgTxt("Usage: 2 inputs and 1 output\n");
}


