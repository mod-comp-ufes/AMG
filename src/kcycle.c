# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include "amg.h"
# include "list.h"
# include "fib_heap.h"
# include "util.h"
# include "solvers.h"
# include "getTime.h"
# include "dpa.h"

int last_i = 0;

void kcycle(int NCL, matrix **A, double **f, double **u, int k, int lmax, double *r, matrix **I_cf, matrix **I_cf_t, double omega, int nr, double tol, int k0, int kw, int i) {

  if(i < (NCL)) {
    
    //if ((i != last_i) || i==0) { // para nao relaxar duas vezes no mesmo nivel no w-cycle
      //SOR_relax(A[i], f[i], u[i], omega, nr);
    //}
    residual(r, f[i], A[i], u[i], 0, NULL, NULL);
    mat_vec(f[i+1], I_cf_t[i], r);
    for (int j=0; j<A[i+1]->n; j++) u[i+1][j] = 0.0;

    if ((((i+1) % k0) !=0 ) || (i == (NCL-1)))
      kcycle(NCL,A,f,u,k,lmax,r,I_cf,I_cf_t,omega,nr,tol,k0,kw,i+1);
    else {
      for(int j = 0; j < kw; j++)
        kcycle(NCL,A,f,u,k,lmax,r,I_cf,I_cf_t,omega,nr,tol,k0,kw,i+1);
    }

    mat_vec_plus(u[i], I_cf[i], u[i+1]);
    //SOR_relax(A[i], f[i], u[i], omega, nr);
  }
  else {
    //SOR_relax(A[i], f[i], u[i], omega, 2*nr);
  }
  GMRES(A[i], f[i], u[i], k, tol, lmax, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, k0, kw);
  last_i = i;

}

void kcycle_precond(int NCL, matrix **A, double **f, double **u, int k, int lmax, double *r, matrix **I_cf, matrix **I_cf_t, double omega, int nr, double tol, int k0, int kw, int i) {
  
  //for (int j=0; j<A[0]->m; j++) u[0][j] = 0.0;

  if(i < (NCL)) {
    
    //if ((i != last_i) || i==0) // para nao relaxar duas vezes no mesmo nivel no w-cycle
      //SOR_relax(A[i], f[i], u[i], omega, nr);
    residual(r, f[i], A[i], u[i], 0, NULL, NULL);
    mat_vec(f[i+1], I_cf_t[i], r);
    for (int j=0; j<A[i+1]->n; j++) u[i+1][j] = 0.0;

    if ((((i+1) % k0) !=0 ) || (i == (NCL-1)))
      kcycle_precond(NCL,A,f,u,k,lmax,r,I_cf,I_cf_t,omega,nr,tol,k0,kw,i+1);
    else {
      for(int j = 0; j < kw; j++)
        kcycle_precond(NCL,A,f,u,k,lmax,r,I_cf,I_cf_t,omega,nr,tol,k0,kw,i+1);
    }

    mat_vec_plus(u[i], I_cf[i], u[i+1]);
    //SOR_relax(A[i], f[i], u[i], omega, nr);
  }
  else {
    //SOR_relax(A[i], f[i], u[i], omega, nr);
    //SOR_relax_rev(A[i], f[i], u[i], omega, nr);
  }
  GMRES(A[i], f[i], u[i], 2, tol, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
  last_i = i;

}

void dpa_kcycle_precond(int NCL, matrix **A, double **f, double **u, int k, int lmax, double *r, int **match, double omega, int nr, double tol, int k0, int kw, int i) {
  
  //for (int j=0; j<A[0]->m; j++) u[0][j] = 0.0;
  //if ((i % 2) == 0)
		//GMRES(A[i], f[i], u[i], 2, tol, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

  if(i < (2*NCL)) {
    
    //if ((i % 2) == 0)
		//if ((i != last_i) || i==0) // para nao relaxar duas vezes no mesmo nivel no w-cycle
			//SOR_relax(A[i], f[i], u[i], omega, nr);
    residual(r, f[i], A[i], u[i], 0, NULL, NULL);
    dpa_fine2coarse(match[i],   r,      A[i]->n  , f[i+1]); // f[i+1] <- r
    dpa_fine2coarse(match[i+1], f[i+1], A[i+1]->n, f[i+2]); // f[i+2] <- f[i+1]
    for (int j=0; j<A[i+1]->n; j++) u[i+1][j] = 0.0;
    for (int j=0; j<A[i+2]->n; j++) u[i+2][j] = 0.0;

    if ((((i+1) % k0) !=0 ) || (i == (2*NCL-1)))
		dpa_kcycle_precond(NCL,A,f,u,k,lmax,r,match,omega,nr,tol,k0,kw,i+2);
    else {
		for(int j = 0; j < kw; j++)
			dpa_kcycle_precond(NCL,A,f,u,k,lmax,r,match,omega,nr,tol,k0,kw,i+2);
    }

    dpa_coarse2fine_pure(match[i+1], u[i+2],   A[i+1]->n, u[i+1]); // u[i-1] <- u[i]
    dpa_coarse2fine     (match[i], u[i+1], A[i]->n, u[i]); // u[i-2] <- u[i-1]
    //if ((i % 2) == 0)
		//SOR_relax(A[i], f[i], u[i], omega, nr);
  }
  else {
    //SOR_relax(A[i], f[i], u[i], omega, nr);
    //SOR_relax_rev(A[i], f[i], u[i], omega, nr);
  }
  if ((i % 2) == 0)
	  GMRES(A[i], f[i], u[i], /*k*/2, tol, /*lmax*/1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
  last_i = i;

}

int AMG_kcycle (matrix *A_o, double *f_o, double *u_o, int krylov, int NCL, double str_thr, int aggress, int refinePass, int truncq, double trunc_fact, double omega, int nr, double tol, int lmax, int k0, int kw) {
	int i, k;
	double **f, **u, *r, delta, t, norm_f;
	matrix **A, **I_cf, **I_cf_t;
	t = AMG_setup(A_o, f_o, u_o, NCL, str_thr, aggress, refinePass, truncq, trunc_fact, &A, &f, &u, &r, &I_cf, &I_cf_t);
	printf("	AMG_kcycle(%d,%.2f,%d,%d,%d,%.2f,%.2f,%d,%e,%d,%d)", NCL, str_thr, aggress, refinePass, truncq, trunc_fact, omega, nr, tol, k0, kw);
	for (i=0; i<NCL+1; i++) printf("	%d	%d", A[i]->m, A[i]->nnz);
	printf("	%.6f", t);
	residual(r, f_o, A_o, u_o, 0, NULL, NULL);
	norm_f = norm_inf(f_o, A_o->m);
	delta = norm_inf(r, A_o->m)/norm_f;
	k = 0;
	t = get_time();
	while (k<lmax && delta>tol) {
		kcycle(NCL, A, f, u, krylov, lmax, r, I_cf, I_cf_t, omega, nr, tol, k0, kw, 0);
		residual(r, f_o, A_o, u_o, 0, NULL, NULL);
		delta = norm_inf(r, A_o->m)/norm_f;
		k++;
		//printf("  %.6f\n", delta);
	}
	t = (get_time() - t)/100.0;
	printf("	%d	%.6f", k, t);
	AMG_destroy(NCL, A, f, u, r, I_cf, I_cf_t);
	return (delta<=tol)? 1: 0;
};

