// #include <R.h>
#include <Rmath.h>

void gcm_nominal(int *y, int *ny, double *x, int *nx, int *bt, int *et, int *lt, double *w, double *r, double *q, double *lambda, double *gamma, double *dist, double *sim, double *ypred)
{
  /*
  y = ny*nt matrix
  x = nx*nt matrix
  w = nx vector
  dist = nt vector
  sim = ny vector
  ypred = ny*nt matrix
  */ 
  // int i,j,t,tt;
  double syp = 0.0;
     
  for(int k=0; k < *lt; k++) {   
	  for(int t=bt[k]; t < et[k]; t++) 
	  {
	    for(int j=0; j < *ny; j++) sim[j] = 0.0;
	    for(int tt=bt[k]-1; tt < t; tt++)
	    {
	      // compute distance
	      dist[tt] = 0.0;
	      for(int i=0; i < *nx; i++)
	      {
	        dist[tt] += w[i] * pow( fabs(x[i + tt * *nx] - x[i + t * *nx ]), *r);
	      }
	      dist[tt] = pow(dist[tt],*q / *r);
	      // compute similarity and add to sim
	      for(int j=0; j < *ny; j++)
	      {
	        if(y[j + tt * *ny] == 1) sim[j] += exp(-1 * *lambda * dist[tt]);
	      }
	    }
	    // make prediction
	    syp = 0.0;
			for(int j=0; j < *ny; j++)
			{
			  ypred[j + t * *ny] = pow(sim[j], *gamma);
			  syp += ypred[j + t * *ny];
			}
	    for(int j=0; j < *ny; j++)
	    {
	      ypred[j + t * *ny] = ypred[j + t * *ny] / syp;
	    }
	  }
  }
}

