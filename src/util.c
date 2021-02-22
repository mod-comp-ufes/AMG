#include <stdlib.h>
#include <math.h>
#include "util.h"
#include "getTime.h"
#include "dpa.h"
#include "ILU/ilup.h"
#include "kcycle.h"

void precond_sol(double *q, matrix *A, double *p, int precond, void * set, SparILU* ILUp) { // solve Mp = q
    precondAMG *amg_s;
    precondDPA *amg_d;
    int i, j, n = A->m;
    double *w = (double *) malloc(n * sizeof (double));
    switch (precond) { // preconditioning
        case 1: // Gauss-Seidel
            for (i = 0; i < n; i++) {
                w[i] = q[i];
                for (j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++) {
                    if (A->col_ind[j] < i) w[i] -= A->val[j] * w[A->col_ind[j]];
                }
                w[i] /= A->diag[i];
            }
            for (i = n - 1; i >= 0; i--) {
                p[i] = w[i];
                for (j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++) {
                    if (A->col_ind[j] > i) p[i] -= A->val[j] * p[A->col_ind[j]];
                }
                p[i] /= A->diag[i];
            }
            break;
        case 2: // AMG
            amg_s = (precondAMG*)set;
            amg_s->f[0] = q;
            amg_s->u[0] = p;
            AMG_Vcycle_precond(amg_s->NCL, amg_s->A, amg_s->f, amg_s->u, amg_s->r, amg_s->I_cf, amg_s->I_cf_t, amg_s->omega, amg_s->nr);
            break;
        case 3: // DPA AMG
            amg_d = (precondDPA*)set;
            amg_d->f[0] = q;
            amg_d->u[0] = p;
            dpa_Vcycle_precond2(amg_d->NCL, amg_d->A, amg_d->f, amg_d->u, amg_d->r, amg_d->omega, amg_d->nr, amg_d->match);
            break;
        case 4: // ILU p 
            ILUp_precond(ILUp, A->n, q, p);
            break;
        case 5: // AMG K ciclo
            amg_s = (precondAMG*)set;
            amg_s->f[0] = q;
            amg_s->u[0] = p;
			kcycle_precond(amg_s->NCL, amg_s->A, amg_s->f, amg_s->u, 20, 1000, amg_s->r, amg_s->I_cf, amg_s->I_cf_t, amg_s->omega, amg_s->nr, 1e-8, amg_s->k0, amg_s->kw, 0);
			break;
        case 6: // DPA AMG K ciclo
            amg_d = (precondDPA*)set;
            amg_d->f[0] = q;
            amg_d->u[0] = p;
			dpa_kcycle_precond(amg_d->NCL, amg_d->A, amg_d->f, amg_d->u, 20, 1000, amg_d->r, amg_d->match, amg_d->omega, amg_d->nr, 1e-8, amg_d->k0, amg_d->kw, 0);
			break;
    }
    free(w);
};

void matvec_product(double *p, matrix *A, double *v, int precond, void * set, SparILU* ILUp) {
    double *q;
    if (precond) {
        q = (double *) malloc(A->m * sizeof (double));
        mat_vec(q, A, v);
        precond_sol(q, A, p, precond, set, ILUp);
        free(q);
    } else {
        mat_vec(p, A, v);
    }
};

void residual(double *r, double *f, matrix *A, double *u, int precond, void * set, SparILU* ILUp) {
    int i, n = A->m;
    matvec_product(r, A, u, precond, set, ILUp);
    for (i = 0; i < n; i++) r[i] = f[i] - r[i];
};

double norm_inf(double *v, int n) {
    int i;
    double max = fabs(v[0]);
    for (i = 1; i < n; i++) {
        if (fabs(v[i]) > max) max = fabs(v[i]);
    }
    return max;
};

double norm_euclid(double *v, int n) {
    int i;
    double sum = 0.0;
    for (i = 0; i < n; i++) sum += v[i] * v[i];
    return sqrt(sum);
};

double ddot(int n, double *dx, double *dy) {
    double ddot = 0.0;
    int i, m;
    if (n <= 0) return ddot;
    m = n % 5;
    i = 0;
    if (m) {
        for (; i < m; i++) ddot += dx[i] * dy[i];
        if (n < 5) return ddot;
    }
    for (; i < n; i += 5) {
        ddot += dx[i] * dy[i] + dx[i + 1] * dy[i + 1] + dx[i + 2] * dy[i + 2] + dx[i + 3] * dy[i + 3] + dx[i + 4] * dy[i + 4];
    }
    return ddot;
};

void daxpy(int n, double da, double *dx, double *dy) {
    int i, m;
    if (n <= 0) return;
    m = n % 4;
    i = 0;
    if (m) {
        for (; i < m; i++) dy[i] += da * dx[i];
    }
    if (n < 4) return;
    for (; i < n; i += 4) {
        dy[i] += da * dx[i];
        dy[i + 1] += da * dx[i + 1];
        dy[i + 2] += da * dx[i + 2];
        dy[i + 3] += da * dx[i + 3];
    }
};

void daypx(int n, double da, double *dx, double *dy) {
    int i, m;
    if (n <= 0) return;
    m = n % 4;
    i = 0;
    if (m) {
        for (; i < m; i++) dy[i] = da * dy[i] + dx[i];
    }
    if (n < 4) return;
    for (; i < n; i += 4) {
        dy[i] = da * dy[i] + dx[i];
        dy[i + 1] = da * dy[i + 1] + dx[i + 1];
        dy[i + 2] = da * dy[i + 2] + dx[i + 2];
        dy[i + 3] = da * dy[i + 3] + dx[i + 3];
    }
};

double dmax(double *v, int n) {
    int i;
    double max = v[0];
    for (i = 1; i < n; i++) {
        if (v[i] > max) max = v[i];
    }
    return max;
};

double dmin(double *v, int n) {
    int i;
    double min = v[0];
    for (i = 1; i < n; i++) {
        if (v[i] < min) min = v[i];
    }
    return min;
};

int compare_arrays(void* vetor1, void* vetor2, int length, int typesize, double tol) {
    int tamanho = length*typesize;
    unsigned char *v1 = (unsigned char*) vetor1;
    unsigned char *v2 = (unsigned char*) vetor2;
    for (int i = 0; i < tamanho; i++)
        if (abs(v1[i] - v2[i]) > tol)
            return 0;
    return 1;
}
