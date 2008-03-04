#include <R.h>
#include <Rmath.h>

/*

NAME: discrBayes
GOAL: computes the online regression weights for a Bayesian logistic regression model

AUTHOR: M. Speekenbrink
VERSION: 0.1
DATE: 06-02-2006
      
ARGUMENTS:
y: vector with y values (0 or 1)
x: as.vector(t(matrix(x))) with x values (0 or 1) (j varies first, then t)
nrx: nrow(x)
ncx: ncol(x)
g: stacked vector of grid values for each dimension
ng: vector of length ncx with number of grid points for each dimension
sw: vector for sums of all combinations of grid points, can be 0 vector of length nsw
nsw: total number of points in multidimensional grid
prior: prior probabilities for each grid point (updates to posterior)
lik: stacked vector for likelihoods, can be 0 vector of length nsw
beta: scale parameter for logistic
weight: stacked vector for output weights, can be 0 vector of length nrx*ncx
pred: vector of lenght ncx for predictions, 
      
RETURNS:
weight: an nx*t matrix with regression weights
pred: a vector of length t with p(y=1)

Algorithm overview:
while t < nt DO
    2. compute likelihood P(y_t|w,x_t)
    3. compute likelihood*prior
    4. normalise prior (is posterior)
    5. compute expected values E[w]
end DO
*/

void discrBayes(int *y, int *x, int *nrx, int *ncx, double *g, int *ng, double *sw, int *nsw, double *prior, double *lik, double *beta, double *weight, double *pred)
{
    int i,j,t,ngt,ngt2;
    int k=0;
    double norm, tweight, tpred, xt;
            
    for(t = 0; t < *nrx; t++) {
        
        for(i = 0; i < *nsw; i++)
              *(sw + i) = 0.0; // initialize sw to 0
        
        ngt = 1; // initialize ngt (counter)
        ngt2 = 0; // initialize ngt2 (counter)
        for(j = 0; j < *ncx; j++) {
            xt = (double) *(x + t * *ncx + j);
            if(xt != 0) {
                for(i = 0; i < *nsw; i++) {
                k = (int) ftrunc(i/ngt) % *(ng + j);
                sw[i] += *(g + k + ngt2) * xt;
                }
            }
            ngt *= *(ng + j);
            ngt2 += *(ng + j);
        }
        norm = 0.0;
        tpred = 0.0;
        for(i = 0; i < *nsw; i++) {
            lik[i] = (*(y + t) * 1/(1+exp(-*beta * *(sw + i)))) + (1-*(y + t)) * (1-1/(1+exp(-*beta * *(sw + i))));
            tpred += *(prior + i) * *(lik + i);
            prior[i] = *(prior + i) * *(lik + i);
            norm += *(prior + i);
        }
        for(i = 0; i < *nsw; i++)
            prior[i] = *(prior + i)/norm;
        
        ngt = 1; // initialize ngt (counter)
        ngt2 = 0;
        for(j = 0; j < *ncx; j++) {
            tweight = 0.0;
            for(i = 0; i < *nsw; i++) {
                  k = (int) ftrunc(i/ngt) % *(ng + j);
                  tweight += *(prior + i) * *(g + ngt2 + k);
            }
            weight[t * *ncx + j] = tweight;
            ngt *= *(ng + j);
            ngt2 += *(ng + j);
        }
  //      tpred = 0.0;
  //      for(i = 0; i < *nsw; i++) {
  //            tpred += *(prior + i) * *(lik + i);
  //      }
        pred[t] = tpred;
    }
}
