#include <R.h>
#include <Rmath.h>

void slfn_logistic_1y(int *y, double *x, int *nx, int *nt, double *eta, double *w)
{
  /*
  y = nt*1 vector
  x = nx*nt matrix
  w = nx*nt matrix
  */ 
 
  int i,t;
  double ypl = 0.0;
  double yp = 0.0;
  
  // make prediction
  for(i=0; i < *nx; i++) {
      ypl += w[i] * x[i];
   }
   yp = 1/(1+exp(-1 * ypl));
    
  for(t=0; t < *nt; t++) {
  	ypl = 0.0;
    for(i=0; i < *nx; i++) {
    
      // update weights
      
      w[i + (t+1) * *nx] = w[i + t * *nx] - eta[i]*(y[t] - yp)*x[i + t * *nx] ;
      // make prediction for next trial
      ypl += w[i + t * *nx] * x[i + t * *nx] ;
    }
    yp = 1/(1+exp(-1 * ypl));
  }

}
