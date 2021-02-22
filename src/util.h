# ifndef util_h
# define util_h
# include "matrix.h"
# include "amg.h"
# include "ILU/ilup.h"

void precond_sol (double *q, matrix *A, double *p, int precond, void * set, SparILU* ILUp);
void matvec_product (double *p, matrix *A, double *v, int precond, void * set, SparILU* ILUp);
void residual (double *r, double *f, matrix *A, double *u, int precond, void * set, SparILU* ILUp);
double norm_inf (double *v, int n);
double norm_euclid (double *v, int n);
double ddot (int n, double *dx, double *dy);
void daxpy (int n, double da, double *dx, double *dy);
void daypx (int n, double da, double *dx, double *dy);
double dmax (double *v, int n);
double dmin (double *v, int n);
int compare_arrays(void* vetor1, void* vetor2, int length, int typesize, double tol);

# endif
