# include "matrix.h"
double error_setup (int NCL, double ***f, double ***u, double *r, int n);
void error_array_sum(double *v1, double *v2, int n);
void error_Vcycle (int NCL, matrix *A, double **f, double **u, double *r, double omega, int nr);
void error_destroy (int NCL, double **f, double **u);
int error(matrix *A_o, double *f_o, double *u_o, int NCL, double omega, int nr, double tol, int lmax);
 