#include <R.h>
#include <Rmath.h>
#include <matrix.h>
#include <R_ext/Applic.h> // for demgg

// the following matrix functions are taken from array.c

// *** BEGIN MATRIX FUNCTIONS
static void matprod(double *x, int nrx, int ncx,
		    double *y, int nry, int ncy, double *z)
{
    char *transa = "N", *transb = "N";
    int i,  j, k;
    double one = 1.0, zero = 0.0;
    LDOUBLE sum;
    Rboolean have_na = FALSE;

    if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
	/* Don't trust the BLAS to handle NA/NaNs correctly: PR#4582
	 * The test is only O(n) here
	 */
	for (i = 0; i < nrx*ncx; i++)
	    if (ISNAN(x[i])) {have_na = TRUE; break;}
	if (!have_na)
	    for (i = 0; i < nry*ncy; i++)
		if (ISNAN(y[i])) {have_na = TRUE; break;}
	if (have_na) {
	    for (i = 0; i < nrx; i++)
		for (k = 0; k < ncy; k++) {
		    sum = 0.0;
		    for (j = 0; j < ncx; j++)
			sum += x[i + j * nrx] * y[j + k * nry];
		    z[i + k * nrx] = sum;
		}
	} else
	    F77_CALL(dgemm)(transa, transb, &nrx, &ncy, &ncx, &one,
			    x, &nrx, y, &nry, &zero, z, &nrx);
    } else /* zero-extent operations should return zeroes */
	for(i = 0; i < nrx*ncy; i++) z[i] = 0;
}

static void symcrossprod(double *x, int nr, int nc, double *z)
{
    char *trans = "T", *uplo = "U";
    double one = 1.0, zero = 0.0;
    int i, j;
    if (nr > 0 && nc > 0) {
        F77_CALL(dsyrk)(uplo, trans, &nc, &nr, &one, x, &nr, &zero, z, &nc);
	for (i = 1; i < nc; i++)
	    for (j = 0; j < i; j++) z[i + nc *j] = z[j + nc * i];
    } else { /* zero-extent operations should return zeroes */
	for(i = 0; i < nc*nc; i++) z[i] = 0;
    }

}

static void crossprod(double *x, int nrx, int ncx,
		      double *y, int nry, int ncy, double *z)
{
    char *transa = "T", *transb = "N";
    double one = 1.0, zero = 0.0;
    if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
        F77_CALL(dgemm)(transa, transb, &ncx, &ncy, &nrx, &one,
			x, &nrx, y, &nry, &zero, z, &ncx);
    } else { /* zero-extent operations should return zeroes */
	int i;
	for(i = 0; i < ncx*ncy; i++) z[i] = 0;
    }
}

// Quick matrix product: no checks!
static void qmatprod(double *x, int nrx, int ncx,
		      double *y, int nry, int ncy, double *z)
{
    char *transa = "N", *transb = "N";
    double one = 1.0, zero = 0.0;
    if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
        F77_CALL(dgemm)(transa, transb, &ncx, &ncy, &nrx, &one,
			x, &nrx, y, &nry, &zero, z, &ncx);
    } else { /* zero-extent operations should return zeroes */
	int i;
	for(i = 0; i < ncx*ncy; i++) z[i] = 0;
    }
}

static void symtcrossprod(double *x, int nr, int nc, double *z)
{
    char *trans = "N", *uplo = "U";
    double one = 1.0, zero = 0.0;
    int i, j;
    if (nr > 0 && nc > 0) {
        F77_CALL(dsyrk)(uplo, trans, &nr, &nc, &one, x, &nr, &zero, z, &nr);
	for (i = 1; i < nr; i++)
	    for (j = 0; j < i; j++) z[i + nr *j] = z[j + nr * i];
    } else { /* zero-extent operations should return zeroes */
	for(i = 0; i < nr*nr; i++) z[i] = 0;
    }

}

static void tcrossprod(double *x, int nrx, int ncx,
		      double *y, int nry, int ncy, double *z)
{
    char *transa = "N", *transb = "T";
    double one = 1.0, zero = 0.0;
    if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
        F77_CALL(dgemm)(transa, transb, &nrx, &nry, &ncx, &one,
			x, &nrx, y, &nry, &zero, z, &nrx);
    } else { /* zero-extent operations should return zeroes */
	int i;
	for(i = 0; i < nrx*nry; i++) z[i] = 0;
    }
}

// *** END MATRIX FUNCTIONS


// function(y,x,A,Q,R,ws,Sigma,nt=nrow(x),nx=ncol(x),lt=1,bt=1,et=nt)

void dlrm_kf_ext(double *y,double *x,double *A, double *Q,double *R, double *ws,
  double *Sigma, int *nt, int *nx, int *ny, int *lt, int *bt, int *et,
  double *w, double *L, double *P, double *H, double *onestep, double *like)
{
  double[] II // II <- diag(nx)
  double[][] t_P, t_L, B, K // t_P <- t_L <- matrix(0.0,nrow=nrow(Q),ncol=ncol(Q))
  //doubleB <- matrix(0.0,ncol=nx)
  //K <- matrix(0.0,ncol=ny,nrow=nx)
  
  for(j=0, j++, j<lt)
  {
    qmatprod(*A,nx,nx,*ws,1,nx,*w[nx*bt[j]]) // w[bt[j],] <- A%*%ws
    // P[bt[j],] <- as.vector(A%*%Sigma%*%t(A) + Q)
    for(i in bt[case]:(et[case]-1))
    {
      B <- matrix(x[i,],ncol=nx)
      t_P <- matrix(P[i,],ncol=nx)
      H[i] <- 1/(B%*%t_P%*%t(B) + R)
      K <- t_P%*%t(B)%*%H[i]
      t_L <- A%*%(II - K%*%B)
      L[i,] <- as.vector(t_L)

      onestep[i] <- B%*%as.matrix(w[i,])
      like <- like + dnorm(y[i,],onestep[i],sqrt(1/H[i]),log=TRUE)

      w[i+1,] <- A%*%K%*%y[i,] + t_L%*%as.matrix(w[i,])
      P[i+1,] <- as.vector(t_L%*%t_P%*%t(A) + Q)
    }
    onestep[et[case]] <- matrix(x[et[case],],ncol=nx)%*%as.matrix(w[et[case],])
  }

  
  
}