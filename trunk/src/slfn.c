// #include <R.h>
#include <Rmath.h>

void slfn(double *y, int *ny, double *x, int *nx, int *nt, double *eta, int *actfun, double *w, double *ypred)
{
  /*
  y = ny*nt matrix
  x = nx*nt matrix
  eta = nx*ny matrix
  actfun: 1 = linear, 2 = logistic, 3 = softmax
  w = nx*ny*nt matrix
  ypred = ny*nt matrix
  */ 
 
  int i,j,t;
  int dim1 = *nx;
  int dim2 = *nx * *ny;
  double ypl[*ny];
  double syp = 0.0;
  
  // make prediction
  for(j=0; j < *ny; j++)
  {
    ypl[j] = 0.0;
    for(i=0; i < *nx; i++) { ypl[j] += w[i + j*dim1] * x[i]; }
    if(*actfun == 1) { ypred[j] = ypl[j]; }
    if(*actfun == 2) { ypred[j] = 1/(1+exp(-1 * ypl[j])); }
    if(*actfun == 3) 
    { 
      ypred[j] = exp(ypl[j]); 
      syp += ypred[j];
    }
  }
  if(*actfun == 3) for(j=0; j < *ny; j++) { ypred[j] = ypred[j]/syp; }
    
  for(t=0; t < *nt; t++) 
  {
    syp = 0.0;
    for(j = 0; j < *ny; j++)
    {
  	  ypl[j] = 0.0;
      for(i=0; i < *nx; i++) 
      {
        // update weights
        w[i + j * dim1 + (t+1) * dim2] = w[i + j * dim1 + t * dim2] - eta[i + j * dim1]*(y[j + t * dim1] - ypred[t * *ny])*x[i + j * dim1 + t * dim2];
        // make prediction for next trial
        ypl[j] += w[i + j * dim1 + t * dim2] * x[i + j * dim1 + t * dim2];
      }
      if(*actfun == 1) { ypred[j + t * *ny] = ypl[j]; }
      if(*actfun == 2) { ypred[j + t * *ny] = 1/(1+exp(-1 * ypl[j])); }
      if(*actfun == 3) 
      { 
        ypred[j + t * *ny] = exp(ypl[j]);
        syp += ypred[j];
      }
    }
    if(*actfun == 3) for(j=0; j < *ny; j++) { ypred[j] = ypred[j]/syp; }
  }
}
