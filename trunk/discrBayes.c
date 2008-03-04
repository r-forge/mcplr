#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stdlib.h>

/*
NAME: discrBayes
GOAL: computes the online regression weights for a Bayesian logistic regression model

AUTHOR: M. Speekenbrink
VERSION: 0.1
DATE: 02-02-2006
      
ARGUMENTS:
nsw: total number of points in multidimensional grid
sw: vector with sums of all combinations of grid points
y: vector with y values (0 or 1)
x: as.vector(t(matrix(x))) with x values (0 or 1) (j varies first, then t)
nrx: nrow(x)
ncx: ncol(x)
g: grid values for each dimension
prior: prior probabilities for each grid point (updates to posterior)

RETURNS:
weight: an nx*t matrix with regression weights

Algorithm overview:
while t < nt DO
    2. compute likelihood P(y_t|w,x_t)
    3. compute likelihood*prior
    4. normalise prior (is posterior)
    5. compute expected values E[w]
end DO
*/

SEXP discrBayes(SEXP y, SEXP x, int *nsw, int *nrx, int *ncx, double *g, double *prior)
{
    int ngt;
    double sw, norm, tweight, *post;
    SEXP weight;
    
    post = (double*) malloc( sizeof(double) * *nsw); // Dynamically assign a vector for post
    assert(post != NULL);
    
    sw = (double*) malloc( sizeof(double) * double(*nsw)); // Dynamically assign a vector for sw
    assert(sw != NULL);
    
    PROTECT(weight = allocMatrix(REALSXP, *nrx, *ncx));
    
    for(t = 0; t < nrx; t++) {
        
        ngt = 1; // initialize ngt (counter)
        for(i = 0; i < nsw; i++)
            sw[i] = 0.0; // initialize sw to 0
        
        for(j = 0; j < ncx; j++) {
            for(i = 0; i < nsw; i++)
                sw[i] += *g[ngt-1 + int(floor(double(i)/double(ngt))) % *ng[j]] * *x[t*nx+j];
        }
        for(i = 0; i < nsw; i++) {
            lik[i] = *y[t] - 1/1+exp(-sw[i]);
            post[i] = *prior[i] * lik[i];
        }
        norm = 0.0;
        for(i = 0; i < nsw; i++)
            norm =+ *prior[i];
        for(i = 0; i < nsw; i++)
            prior[i] = *prior[i]/norm;

        tweight = 0.0;
        for(j = 0; j < ncx; j++) {
            for(i = 0; i < nsw; i++)
                tweight += *prior[i] * *g[ngt-1 + floor(i/ngt) % *ng[j]];
            weight[t*nx+j] = tweight;
        }
    }
    UNPROTECT(1);
    free(post);
    return(weight)
}
