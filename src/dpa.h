# ifndef dpa_h
# define dpa_h

#include "matrix.h"

typedef struct {
	int NCL, nr, **match;
	double **f, **u, *r, omega;
	matrix **A;
	int k0, kw;
} precondDPA;

double dpa_setup(matrix *A_o, double *f_o, double *u_o, int NCL, matrix ***A, double ***f, double ***u, double **r, int ***match, double beta);
int		dpa_AMG (matrix *A_o, double *f_o, double *u_o, int NCL, double omega, int nr, double tol, int lmax, double beta);
void	dpa_Vcycle_precond(int NCL, matrix **A, double **f, double **u, double *r, double omega, int nr, int **match);
void	dpa_Vcycle_precond2(int NCL, matrix **A, double **f, double **u, double *r, double omega, int nr, int **match);
void	dpa_destroy(int NCL, matrix **A, double **f, double **u, double *r, int **match);
void	dpa_fine2coarse(int * match, double * fine, int n, double * coarse);
void	dpa_coarse2fine(int * match, double * coarse, int n, double * fine);
void	dpa_coarse2fine_pure(int * match, double * coarse, int n, double * fine);

# endif
