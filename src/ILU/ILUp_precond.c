#include "ilup.h"

int ILUp_precond(SparILU *ILUp, int neq, double *p, double *z)
{
	/*--------------------------------------------------------------------
	 performs a forward followed by a backward solve
	  for LU matrix as produced by ilup
	  y  = right-hand-side
	  x  = solution on return
	  lu = LU matrix as produced by ilup.
	------------------------------------------------------------------*/
	int i, j, nzcount, *ja;
	double *D, *ma, *x, *y;
	SparMAT *L, *U;
	SparILU *lu;	

	lu = ILUp;
	L = lu->L;
	U = lu->U;
	D = lu->D;
	x = z;	
	y = p;

	/* Block L solve */
	for( i = 0; i < neq; i++ ) {
		x[i] = y[i];
		nzcount = L->nzcount[i];
		ja = L->ja[i];
		ma = L->ma[i];
		for( j = 0; j < nzcount; j++ ) {
			x[i] -= x[ja[j]] * ma[j];
		}
	}
	/* Block -- U solve */
	for( i = neq-1; i >= 0; i-- ) {
		nzcount = U->nzcount[i];
		ja = U->ja[i];
		ma = U->ma[i];
		for( j = 0; j < nzcount; j++ ) {
			x[i] -= x[ja[j]] * ma[j];
		}
		x[i] *= D[i];
	}

	return 0;
}

