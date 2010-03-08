// #include <R.h>
#include <Rmath.h>

void gcm_nominal(int *y, int *ny, double *x, int *nx, int *nt, double *w, double *r, double *q, double *lambda, double *gamma, double *dist, double *sim, double *ypred)
{
  /*
  y = ny*nt matrix
  x = nx*nt matrix
  w = nx vector
  dist = nt vector
  sim = ny vector
  ypred = ny*nt matrix
  */ 
  int i,j,t,tt;
  double syp = 0.0;
     
  for(t=1; t < *nt; t++) 
  {
    for(j=0; j < *ny; j++) sim[j] = 0.0;
    for(tt=0; tt < t; tt++)
    {
      // compute distance
      dist[tt - 1] = 0.0;
      for(i=0; i < *nx; i++) dist[tt-1] += w[i] * pow( fabs(x[i + tt * *nx] - x[i + t * *nx ]), *r);
      dist[tt - 1] = pow(dist[tt-1],1 / *r);
      // compute similarity and add to sim
      for(j=0; j < *ny; j++)
      {
        if(y[j + tt * *ny] == 1) sim[j] += exp( - *lambda * pow(dist[tt - 1],*q));
      }
    }
    // make prediction
    syp = 0.0;
	for(j=0; j < *ny; j++)
	{
	  ypred[j + t * *ny] = pow(sim[j], *gamma);
	  syp += ypred[j + t * *ny];
	}
    for(j=0; j < *ny; j++) ypred[j + t * *ny] = ypred[j + t * *nt]/syp;
  }
}
