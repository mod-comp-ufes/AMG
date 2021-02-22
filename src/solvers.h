# ifndef solvers_h
# define solvers_h
# include "matrix.h"

void SOR_relax (matrix *A, double *f, double *u, double omega, int iter);
void SOR_relax_rev (matrix *A, double *f, double *u, double omega, int iter);
int SOR (matrix *A, double *f, double *u, double omega, double tol, int lmax);
int CG (matrix *A, double *f, double *u, double tol, int lmax);
int GMRES (matrix *A, double *f, double *u, int k, double tol, int lmax, int precond, int NCL, double str_thr, int aggress, int refinePass, int truncq, double trunc_fact, double omega, int nr, int ILU_p, int k0, int kw);
int LCD (matrix *A, double *f, double *u, int k, double tol, int lmax, int precond, int NCL, double str_thr, int aggress, int refinePass, int truncq, double trunc_fact, double omega, int nr, int ILU_p);

# endif
