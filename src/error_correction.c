# include <stdlib.h>
# include <stdio.h>
# include <float.h>
# include <math.h>
# include "getTime.h"
# include "matrix.h"
# include "solvers.h"
# include "util.h"

double error_setup (int NCL, double ***f, double ***u, double **r, int n) {
	// SETUP PHASE
	double t = get_time();
	
	*f = (double **) malloc((NCL+1)*sizeof(double *));
	*u = (double **) malloc((NCL+1)*sizeof(double *));

	for (int i=1; i<NCL+1; i++) {
		(*f)[i] = (double *) malloc(n*sizeof(double));
		(*u)[i] = (double *) malloc(n*sizeof(double));
	}
	t = (get_time() - t)/100.0;
	*r = (double *) malloc(n*sizeof(double));
	return t;
};

void error_array_sum(double *v1, double *v2, int n) {
    for(int i = 0; i < n; i++) {
        v1[i] += v2[i];
    }
}

void error_Vcycle (int NCL, matrix *A, double **f, double **u, double *r, double omega, int nr) {

	int i, j;
	for (i=0; i<NCL; i++) {
		SOR_relax(A, f[i], u[i], omega, nr);
		residual(r, f[i], A, u[i], 0, 0, 0);
		f[i+1] = r;
		for (j=0; j<A->n; j++) u[i+1][j] = 0.0;
	}
	SOR_relax(A, f[i], u[i], omega, 2*nr);
	for (i--; i>=0; i--) {
 		// u[i] <- u[i+1]
                error_array_sum(u[i],u[i+1],A->n);
		SOR_relax(A, f[i], u[i], omega, nr);
	}
};

void error_destroy (int NCL, double **f, double **u) {
	for (int i=1; i<NCL+1; i++) { free(f[i]); free(u[i]); }
	free(f); free(u);
};

int error(matrix *A_o, double *f_o, double *u_o, int NCL, double omega, int nr, double tol, int lmax) {
	int k;
	double **f, **u, *r, delta, t, norm_f;

	t = error_setup (NCL, &f, &u, &r, A_o->n);
        f[0] = f_o;
        u[0] = u_o;
        
	printf("\n\tsetup t: %f",t);
        
        //printf("coe\n"); exit(0);
	residual(r, f_o, A_o, u_o, 0, 0, 0);
	norm_f = norm_inf(f_o, A_o->m);
	delta = norm_inf(r, A_o->m)/norm_f;
	k = 0;
	t = get_time();
	while (k<lmax && delta>tol) {
		error_Vcycle(NCL, A_o, f, u, r, omega, nr);
		residual(r, f_o, A_o, u_o, 0, 0, 0); 
		delta = norm_inf(r, A_o->m)/norm_f;
		k++;
		//printf("  %.6f\n", delta);
	}
	t = (get_time() - t)/100.0;
	printf("\n\titer: %d \texec t: %.6f", k, t);
	return (delta<=tol)? 1: 0;
};
