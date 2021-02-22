# include "matrix.h"
void kcycle         	(int NCL, matrix **A, double **f, double **u, int k, int lmax, double *r, matrix **I_cf, matrix **I_cf_t, double omega, int nr, double tol, int k0, int kw, int i);
void kcycle_precond		(int NCL, matrix **A, double **f, double **u, int k, int lmax, double *r, matrix **I_cf, matrix **I_cf_t, double omega, int nr, double tol, int k0, int kw, int i);
void dpa_kcycle_precond	(int NCL, matrix **A, double **f, double **u, int k, int lmax, double *r, int **match, double omega, int nr, double tol, int k0, int kw, int i);
int AMG_kcycle     		(matrix *A_o, double *f_o, double *u_o, int krylov, int NCL, double str_thr, int aggress, int refinePass, int truncq, double trunc_fact, double omega, int nr, double tol, int lmax, int k0, int kw);
