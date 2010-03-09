// #include <R.h>
#include <Rmath.h>

void slfn(double *y, int *ny, double *x, int *nx, int *bt, int *et, int *lt, double *eta, int *actfun, double *w, double *ypred)
{
  /*
  y = ny*nt matrix
  x = nx*nt matrix
  eta = nx*ny matrix
  actfun: 1 = linear, 2 = logistic, 3 = softmax
  w = nx*ny*nt matrix
  ypred = ny*nt matrix
  */ 
 
  //int i,j,t;
  int dim1 = *nx;
  int dim2 = *nx * *ny;
  double ypl[*ny];
  double syp = 0.0;
  // double err = 0.0;
  for(int k=0; k < *lt; k++) 
  {
	  for(int t=(bt[k]-1); t < et[k]; t++) 
	  {
	    // make prediction
	    syp = 0.0;
	    for(int j=0; j < *ny; j++)
	    {
			  ypl[j] = 0.0;
			  for(int i=0; i < *nx; i++)
			  { 
			  	ypl[j] += w[i + j * dim1 + t * dim2] * x[i + t * *nx];
			  }
			  if(*actfun == 1) ypred[j + t * *ny] = ypl[j];
			  if(*actfun == 2) ypred[j + t * *ny] = 1 / (1+exp(- ypl[j]));
			  if(*actfun == 3) 
			  { 
		    	ypred[j + t * *ny] = exp(ypl[j]); 
		      syp += exp(ypl[j]);
			  }
			}
	    if(*actfun == 3) 
	    {
	    	for(int j=0; j < *ny; j++) 
	    	{
	    	  ypred[j + t * *ny] = ypred[j + t * *ny]/syp;
	    	}
	    }
	    // update weights
	    if(t + 1 < et[k])
	    {
		    for(int j = 0; j < *ny; j++)
		    {
		      for(int i = 0; i < *nx; i++) 
		      {
		        w[i + j * dim1 + (t+1) * dim2] = w[i + j * dim1 + t * dim2] + eta[i + j * dim1] * (y[j + t * *ny] - ypred[j + t * *ny]) * x[i + t * *nx];
		      }
		    }
	    }
	  }
  }
}

// This function includes momentum through parameter beta

void slfn_b(double *y, int *ny, double *x, int *nx, int *bt, int *et, int *lt, double *eta, double *beta, int *actfun, double *w, double *ypred)
{
  /*
  y = ny*nt matrix
  x = nx*nt matrix
  eta = nx*ny matrix
  actfun: 1 = linear, 2 = logistic, 3 = softmax
  w = nx*ny*nt matrix
  ypred = ny*nt matrix
  */ 
 
  //int i,j,t;
  int dim1 = *nx;
  int dim2 = *nx * *ny;
  double ypl[*ny];
  double syp = 0.0;
  // double err = 0.0;
  for(int k=0; k < *lt; k++) 
  {
	  for(int t=(bt[k]-1); t < et[k]; t++) 
	  {
	    // make prediction
	    syp = 0.0;
	    for(int j=0; j < *ny; j++)
	    {
			  ypl[j] = 0.0;
			  for(int i=0; i < *nx; i++)
			  { 
			  	ypl[j] += w[i + j * dim1 + t * dim2] * x[i + t * *nx];
			  }
			  if(*actfun == 1) ypred[j + t * *ny] = ypl[j];
			  if(*actfun == 2) ypred[j + t * *ny] = 1 / (1+exp(- ypl[j]));
			  if(*actfun == 3) 
			  { 
		    	ypred[j + t * *ny] = exp(ypl[j]); 
		      syp += exp(ypl[j]);
			  }
			}
	    if(*actfun == 3) 
	    {
	    	for(int j=0; j < *ny; j++) 
	    	{
	    	  ypred[j + t * *ny] = ypred[j + t * *ny]/syp;
	    	}
	    }
	    // update weights
	    if(t + 1 < et[k])
	    {
		    for(int j = 0; j < *ny; j++)
		    {
		      for(int i = 0; i < *nx; i++) 
		      {
		        w[i + j * dim1 + (t+1) * dim2] = w[i + j * dim1 + t * dim2] + eta[i + j * dim1] * (y[j + t * *ny] - ypred[j + t * *ny]) * x[i + t * *nx];
		        if(t - bt[k] > -1)
		        {
		          w[i + j * dim1 + (t+1) * dim2] += beta[i + j * dim1] * (w[i + j * dim1 + t * dim2] - w[i + j * dim1 + (t-1) * dim2]);
		        }
		      }
		    }
	    }
	  }
  }
}


// This function includes decreasing eta as eta = eta/t^alpha

void slfn_a(double *y, int *ny, double *x, int *nx, int *bt, int *et, int *lt, double *eta, double *alpha, int *actfun, double *w, double *ypred)
{
  /*
  y = ny*nt matrix
  x = nx*nt matrix
  eta = nx*ny matrix
  actfun: 1 = linear, 2 = logistic, 3 = softmax
  w = nx*ny*nt matrix
  ypred = ny*nt matrix
  */ 
 
  //int i,j,t;
  int dim1 = *nx;
  int dim2 = *nx * *ny;
  double ypl[*ny];
  double syp = 0.0;
  // double err = 0.0;
  for(int k=0; k < *lt; k++) 
  {
	  for(int t=(bt[k]-1); t < et[k]; t++) 
	  {
	    // make prediction
	    syp = 0.0;
	    for(int j=0; j < *ny; j++)
	    {
			  ypl[j] = 0.0;
			  for(int i=0; i < *nx; i++)
			  { 
			  	ypl[j] += w[i + j * dim1 + t * dim2] * x[i + t * *nx];
			  }
			  if(*actfun == 1) ypred[j + t * *ny] = ypl[j];
			  if(*actfun == 2) ypred[j + t * *ny] = 1 / (1+exp(- ypl[j]));
			  if(*actfun == 3) 
			  { 
		    	ypred[j + t * *ny] = exp(ypl[j]); 
		      syp += exp(ypl[j]);
			  }
			}
	    if(*actfun == 3) 
	    {
	    	for(int j=0; j < *ny; j++) 
	    	{
	    	  ypred[j + t * *ny] = ypred[j + t * *ny]/syp;
	    	}
	    }
	    // update weights
	    if(t + 1 < et[k])
	    {
		    for(int j = 0; j < *ny; j++)
		    {
		      for(int i = 0; i < *nx; i++) 
		      {
		        w[i + j * dim1 + (t+1) * dim2] = w[i + j * dim1 + t * dim2] +  (eta[i + j * dim1] * pow(t - bt[k] + 2, -*alpha)) * (y[j + t * *ny] - ypred[j + t * *ny]) * x[i + t * *nx];
		      }
		    }
	    }
	  }
  }
}


// this function includes decreasing learning rate and momentum
void slfn_ab(double *y, int *ny, double *x, int *nx, int *bt, int *et, int *lt, double *eta, double *alpha, double *beta, int *actfun, double *w, double *ypred)
{
  /*
  y = ny*nt matrix
  x = nx*nt matrix
  eta = nx*ny matrix
  actfun: 1 = linear, 2 = logistic, 3 = softmax
  w = nx*ny*nt matrix
  ypred = ny*nt matrix
  */ 
 
  //int i,j,t;
  int dim1 = *nx;
  int dim2 = *nx * *ny;
  double ypl[*ny];
  double syp = 0.0;
  // double err = 0.0;
  for(int k=0; k < *lt; k++) 
  {
	  for(int t=(bt[k]-1); t < et[k]; t++) 
	  {
	    // make prediction
	    syp = 0.0;
	    for(int j=0; j < *ny; j++)
	    {
			  ypl[j] = 0.0;
			  for(int i=0; i < *nx; i++)
			  { 
			  	ypl[j] += w[i + j * dim1 + t * dim2] * x[i + t * *nx];
			  }
			  if(*actfun == 1) ypred[j + t * *ny] = ypl[j];
			  if(*actfun == 2) ypred[j + t * *ny] = 1 / (1+exp(- ypl[j]));
			  if(*actfun == 3) 
			  { 
		    	ypred[j + t * *ny] = exp(ypl[j]); 
		      syp += exp(ypl[j]);
			  }
			}
	    if(*actfun == 3) 
	    {
	    	for(int j=0; j < *ny; j++) 
	    	{
	    	  ypred[j + t * *ny] = ypred[j + t * *ny]/syp;
	    	}
	    }
	    // update weights
	    if(t + 1 < et[k])
	    {
		    for(int j = 0; j < *ny; j++)
		    {
		      for(int i = 0; i < *nx; i++) 
		      {
		        w[i + j * dim1 + (t+1) * dim2] = w[i + j * dim1 + t * dim2] +  (eta[i + j * dim1] * pow(t - bt[k] + 2, -*alpha)) * (y[j + t * *ny] - ypred[j + t * *ny]) * x[i + t * *nx];
		        if(t - bt[k] > -1)
		        {
		          w[i + j * dim1 + (t+1) * dim2] += beta[i + j * dim1] * (w[i + j * dim1 + t * dim2] - w[i + j * dim1 + (t-1) * dim2]);
		        }
		      }
		    }
	    }
	  }
  }
}

