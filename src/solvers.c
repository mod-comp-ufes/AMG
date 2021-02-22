#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "solvers.h"
#include "util.h"
#include "amg.h"
#include "getTime.h"
#include "dpa.h"
#include "ILU/ilup.h"

void SOR_relax(matrix *A, double *f, double *u, double omega, int iter) {
    int i, j, k, n = A->m;
    double sum;
    for (k = 0; k < iter; k++) {
        for (i = 0; i < n; i++) {
            sum = 0.0;
            for (j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++) {
                if (A->col_ind[j] != i) {
                    sum += A->val[j] * u[A->col_ind[j]];
                }
            }
            u[i] += omega * ((f[i] - sum) / A->diag[i] - u[i]);
        }
    }
};

void SOR_relax_rev(matrix *A, double *f, double *u, double omega, int iter) {
    int i, j, k, n = A->m;
    double sum;
    for (k = 0; k < iter; k++) {
        for (i = n - 1; i >= 0; i--) {
            sum = 0.0;
            for (j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++) {
                if (A->col_ind[j] != i) {
                    sum += A->val[j] * u[A->col_ind[j]];
                }
            }
            u[i] += omega * ((f[i] - sum) / A->diag[i] - u[i]);
        }
    }
};

int SOR(matrix *A, double *f, double *u, double omega, double tol, int lmax) {
    int k, n = A->m;
    double *r, delta, t, norm_f;
    r = (double *) malloc(n * sizeof (double));
    //norm_f = norm_inf(f, n);
    residual(r, f, A, u, 0, 0, NULL);
    delta = norm_inf(r, n); //norm_f;
    t = get_time();
    for (k = 0; k < lmax && delta > tol; k++) {
        SOR_relax(A, f, u, omega, 1);
        residual(r, f, A, u, 0, 0, NULL);
        delta = norm_inf(r, n); //norm_f;
        //printf("  %.6f\n", delta);
    }
    t = (get_time() - t) / 100.0;
    printf("	SOR(%.2f,%e,%d)	%d	%d	%d	%.6f", omega, tol, lmax, n, A->nnz, k, t);
    free(r);
    return (delta <= tol) ? 1 : 0;
};

int CG(matrix *A, double *f, double *u, double tol, int lmax) {
    int k, n = A->m;
    double *v, *r, *p, delta, delta_pr, alpha, beta, t, tol2;
    v = (double *) malloc(n * sizeof (double));
    r = (double *) malloc(n * sizeof (double));
    p = (double *) malloc(n * sizeof (double));
    residual(r, f, A, u, 0, 0, NULL);
    for (k = 0; k < n; k++) v[k] = r[k];
    delta = ddot(n, r, r);
    tol2 = tol * tol*delta;
    t = get_time();
    for (k = 0; k < lmax && delta > tol2; k++) {
        mat_vec(p, A, v);
        alpha = delta / ddot(n, v, p);
        daxpy(n, alpha, v, u);
        daxpy(n, -alpha, p, r);
        delta_pr = delta;
        delta = ddot(n, r, r);
        beta = delta / delta_pr;
        daypx(n, beta, r, v);
        //printf("  %.6f %.6f\n", delta, tol2);
    }
    t = (get_time() - t) / 100.0;
    printf("	CG(%e,%d)	%d	%d	%d	%.6f", tol, lmax, n, A->nnz, k, t);
    free(v);
    free(r);
    free(p);
    return (delta <= tol2) ? 1 : 0;
};

int GMRES(matrix *A, double *f, double *u, int k, double tol, int lmax, int precond, int NCL, double str_thr, int aggress, int refinePass, int truncq, double trunc_fact, double omega, int nr, int ILU_p, int k0, int kw) {
    int ind, i, j, l, n = A->m, *p;
    SparILU * ILUp = (SparILU*)malloc(sizeof(SparILU));
    SparMAT * sparMAT = (SparMAT*)malloc(sizeof(SparMAT));
    MAT * Ailu = (MAT*)malloc(sizeof(MAT));
    double **U, **H, *e, *y, *c, *s, *f2, tol2, delta, aux, r, t;
    void*amg_s;
    U = (double **) malloc((k + 1) * sizeof (double *));
    for (ind = 0; ind < k + 1; ind++) U[ind] = (double *) calloc(n, sizeof (double));
    H = (double **) malloc((k + 1) * sizeof (double *));
    for (ind = 0; ind < k + 1; ind++) H[ind] = (double *) malloc(k * sizeof (double));
    e = (double *) malloc((k + 1) * sizeof (double));
    y = (double *) malloc(k * sizeof (double));
    c = (double *) malloc(k * sizeof (double));
    s = (double *) malloc(k * sizeof (double));
    if ((precond == 2) || (precond == 5)) {
        t = get_time();
        amg_s = malloc(sizeof (precondAMG));
        ((precondAMG *) amg_s)->NCL = NCL;
        ((precondAMG *) amg_s)->omega = omega;
        ((precondAMG *) amg_s)->nr = nr;
        ((precondAMG *) amg_s)->k0 = k0;
        ((precondAMG *) amg_s)->kw = kw;
        AMG_setup(A, f, u, NCL, str_thr, aggress, refinePass, truncq, trunc_fact, &(((precondAMG *) amg_s)->A), &(((precondAMG *) amg_s)->f), &(((precondAMG *) amg_s)->u), &(((precondAMG *) amg_s)->r), &(((precondAMG *) amg_s)->I_cf), &(((precondAMG *) amg_s)->I_cf_t));
    } else if ((precond == 3) || (precond == 6)) {
        t = get_time();
        amg_s = malloc(sizeof (precondDPA));
        ((precondDPA *) amg_s)->NCL = NCL;
        ((precondDPA *) amg_s)->omega = omega;
        ((precondDPA *) amg_s)->nr = nr;
        ((precondDPA *) amg_s)->k0 = k0;
        ((precondDPA *) amg_s)->kw = kw;
        dpa_setup(A, f, u, NCL, &(((precondDPA *) amg_s)->A), &(((precondDPA *) amg_s)->f), &(((precondDPA *) amg_s)->u), &(((precondDPA *) amg_s)->r), &(((precondDPA *) amg_s)->match), str_thr);
    } else if (precond == 4) {
        t = get_time();
        // ILU setup
        ILUp_precond_setup(A, &ILUp, &sparMAT, &Ailu, f, ILU_p);
    }
    t = get_time();
    if (precond && (precond != 4)) { // modify right-hand side
        f2 = (double *) calloc(n, sizeof (double));
        precond_sol(f, A, f2, precond, amg_s, ILUp);
        for (i = 0; i < n; i++)
            f[i] = f2[i];
        free(f2);
    }
    tol2 = tol * norm_euclid(f, n);
    l = 0;
    do {
        i = 0;
        residual(U[i], f, A, u, precond, amg_s, ILUp);
        e[i] = norm_euclid(U[i], n);
        for (aux = 1.0 / e[i], ind = 0; ind < n; ind++) U[i][ind] *= aux;
        delta = e[i];

        for (i; i < k && delta > tol2; i++, l++) {

            matvec_product(U[i + 1], A, U[i], precond, amg_s, ILUp);
            // Gram-Schmidt orthogonalization
            for (j = 0; j <= i; j++) {
                H[j][i] = ddot(n, U[i + 1], U[j]);
                daxpy(n, -H[j][i], U[j], U[i + 1]);
            }
            H[i + 1][i] = norm_euclid(U[i + 1], n);
            for (aux = 1.0 / H[i + 1][i], ind = 0; ind < n; ind++) U[i + 1][ind] *= aux;
            // QR algorithm
            for (j = 0; j <= i - 1; j++) {
                aux = H[j][i];
                H[j][i] = c[j] * H[j][i] + s[j] * H[j + 1][i];
                H[j + 1][i] = -s[j] * aux + c[j] * H[j + 1][i];
            }
            r = sqrt(H[i][i] * H[i][i] + H[i + 1][i] * H[i + 1][i]);
            c[i] = H[i][i] / r;
            s[i] = H[i + 1][i] / r;
            H[i][i] = r; // H[i+1][i] = 0.0;
            e[i + 1] = -s[i] * e[i];
            e[i] = c[i] * e[i];
            delta = fabs(e[i + 1]);
        }
        i--;
        for (j = i; j >= 0; j--) {
            for (aux = 0.0, ind = j + 1; ind <= i; ind++) aux += H[j][ind] * y[ind];
            y[j] = (e[j] - aux) / H[j][j];
        }
        for (ind = 0; ind < n; ind++) {
            for (j = 0; j <= i; j++) u[ind] += y[j] * U[j][ind];
        }
    } while (l < lmax && delta > tol2);

	if (lmax > 2) {
		printf("\tGMRES(%d,%e,%d,%d)	%d	%d\n", k, tol, lmax, precond, n, A->nnz);
		printf("\nAMG (%d, %.1f, %d, %.e, %d, %.2f)", NCL, omega, nr, tol, lmax, str_thr);
	}
	printf("\titer: %d \t\texec: %.6f ", l, t);
	if ((precond == 2) || (precond == 5)) {
		AMG_destroy(((precondAMG *) amg_s)->NCL, ((precondAMG *) amg_s)->A, ((precondAMG *) amg_s)->f, ((precondAMG *) amg_s)->u, ((precondAMG *) amg_s)->r, ((precondAMG *) amg_s)->I_cf, ((precondAMG *) amg_s)->I_cf_t);
		free(((precondAMG *) amg_s));
	} else if ((precond == 3) || (precond == 6)) {
		dpa_destroy(((precondDPA *) amg_s)->NCL, ((precondDPA *) amg_s)->A, ((precondDPA *) amg_s)->f, ((precondDPA *) amg_s)->u, ((precondDPA *) amg_s)->r, ((precondDPA *) amg_s)->match);
		free(((precondDPA *) amg_s));
	} else if (precond == 4) {
		free(ILUp);
	}
    free(s);
    free(c);
    free(y);
    free(e);
    for (ind = 0; ind < k + 1; ind++) free(H[ind]);
    free(H);
    for (ind = 0; ind < k + 1; ind++) free(U[ind]);
    free(U);
    return (delta <= tol2) ? 1 : 0;
};

int LCD(matrix *A, double *f, double *u, int k, double tol, int lmax, int precond, int NCL, double str_thr, int aggress, int refinePass, int truncq, double trunc_fact, double omega, int nr, int ILU_p) {
    int ind, i, j, l, n = A->m, *p;
    double **P, **Q, *PQ, *r, tol2, delta, alfa, beta, t;
    precondAMG *amg_s;
    SparILU *ILUp;
    P = (double **) malloc((k + 1) * sizeof (double *));
    for (ind = 0; ind < k + 1; ind++) P[ind] = (double *) malloc(n * sizeof (double));
    Q = (double **) malloc((k + 1) * sizeof (double *));
    for (ind = 0; ind < k + 1; ind++) Q[ind] = (double *) malloc(n * sizeof (double));
    PQ = (double *) malloc(k * sizeof (double));
    r = (double *) malloc(n * sizeof (double));
    t = get_time();
    if (precond == 2) {
        amg_s = (precondAMG *) malloc(sizeof (precondAMG));
        amg_s->NCL = NCL;
        amg_s->omega = omega;
        amg_s->nr = nr;
        AMG_setup(A, f, u, amg_s->NCL, str_thr, aggress, refinePass, truncq, trunc_fact, &(amg_s->A), &(amg_s->f), &(amg_s->u), &(amg_s->r), &(amg_s->I_cf), &(amg_s->I_cf_t));
    }
    if (precond) {
        precond_sol(f, A, r, precond, amg_s, ILUp);
        for (i = 0; i < n; i++) f[i] = r[i];
    } // modify right-hand side
    residual(r, f, A, u, precond, amg_s, ILUp);
    delta = norm_euclid(r, n);
    for (ind = 0; ind < n; ind++) P[0][ind] = r[ind];
    tol2 = tol * norm_euclid(f, n);
    l = 0;
    do {
        i = 0;
        matvec_product(Q[i], A, P[i], precond, amg_s, ILUp);
        for (i; i < k && delta > tol2; i++, l++) {
            PQ[i] = ddot(n, P[i], Q[i]);
            alfa = ddot(n, P[i], r) / PQ[i];
            daxpy(n, alfa, P[i], u);
            daxpy(n, -alfa, Q[i], r);
            for (ind = 0; ind < n; ind++) P[i + 1][ind] = r[ind];
            matvec_product(Q[i + 1], A, P[i + 1], precond, amg_s, ILUp);
            for (j = 0; j <= i; j++) {
                beta = -ddot(n, P[j], Q[i + 1]) / PQ[j];
                daxpy(n, beta, P[j], P[i + 1]);
                daxpy(n, beta, Q[j], Q[i + 1]);
            }
            delta = norm_euclid(r, n);
        }
        for (ind = 0; ind < n; ind++) P[0][ind] = P[i][ind];
    } while (l < lmax && delta > tol2);
    printf("	LCD(%d,%e,%d,%d)	%d	%d	%d	%.6f", k, tol, lmax, precond, n, A->nnz, l, t);
    printf("\nAMG (%d, %.1f, %d, %.e, %d, %.2f)", NCL, omega, nr, tol, lmax, str_thr);
    if (precond == 2) {
        AMG_destroy(amg_s->NCL, amg_s->A, amg_s->f, amg_s->u, amg_s->r, amg_s->I_cf, amg_s->I_cf_t);
        free(amg_s);
    } else if (precond == 3) {
        dpa_destroy(((precondDPA *) amg_s)->NCL, ((precondDPA *) amg_s)->A, ((precondDPA *) amg_s)->f, ((precondDPA *) amg_s)->u, ((precondDPA *) amg_s)->r, ((precondDPA *) amg_s)->match);
        free(((precondDPA *) amg_s));
    }
    else if (precond == 4) {
        free(ILUp);
    }
    free(r);
    free(PQ);
    for (ind = 0; ind < k + 1; ind++) free(Q[ind]);
    free(Q);
    for (ind = 0; ind < k + 1; ind++) free(P[ind]);
    free(P);
    return (delta <= tol2) ? 1 : 0;
};
