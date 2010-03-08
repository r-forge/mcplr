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
     
  for(t=0; t < *nt; t++) 
  {
    // make prediction
    syp = 0.0;
    for(j=0; j < *ny; j++)
    {
	  ypl[j] = 0.0;
	  for(i=0; i < *nx; i++) ypl[j] += w[i + j * dim1 + t * dim2] * x[i + t * *nx];
	  if(*actfun == 1) ypred[j + t * *ny] = ypl[j];
	  if(*actfun == 2) ypred[j + t * *ny] = 1/(1+exp(-1 * ypl[j]));
	  if(*actfun == 3) 
	  { 
    	ypred[j + t * *ny] = exp(ypl[j]); 
        syp += exp(ypl[j]);
	  }
	}
    if(*actfun == 3) for(j=0; j < *ny; j++) ypred[j + t * *ny] = ypred[j + t * *ny]/syp;
  
    // update weights
    for(j = 0; j < *ny; j++)
    {
      for(i=0; i < *nx; i++) 
      {
        w[i + j * dim1 + (t+1) * dim2] = w[i + j * dim1 + t * dim2] + eta[i + j * dim1]*(y[j + t * *ny] - ypred[j + t * *ny])*x[i + t * *nx];
      }
    }
  }
}
