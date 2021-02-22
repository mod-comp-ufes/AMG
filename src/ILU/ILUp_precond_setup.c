#include "ilup.h"
#include "../matrix.h"

void ILUp(SparMAT* csmat, SparILU* lu, int p);

int ILUp_precond_setup(matrix * matrix, SparILU ** sparILU, SparMAT ** sparMAT, MAT ** Ailu, double *F, int p) {
    /*level of fill-in */
    MAT *A;
    SparMAT *mat;
    SparILU *lu;

    int ierr;

    lu = calloc(1, sizeof (SparILU));
    mat = calloc(1, sizeof (SparMAT));
    A = calloc(1, sizeof (MAT));

    A->n  = matrix->n;
    A->nz = matrix->nnz;
    A->AA = matrix->val;
    A->JA = matrix->col_ind;
    A->IA = matrix->row_ptr;
    A->D  = matrix->diag;

    /*Convert struct CSR in SPARMAT */
    CSRto_SPARMAT_setup(A, mat);

    SPARILU_setup(lu, A->n);

    /* symbolic factorization to calculate level of fill index arrays */
    if ((ierr = LEVEL_OF_FILL_lu(mat, lu, p)) != 0) {
        printf("Error: LEVEL OF FILL\n");
        exit(1);
    }
    (*sparILU)    = lu;
    (*sparMAT) = mat;
    (*Ailu)    = A;

    ILUp(mat, lu, p);

    (*sparILU) = lu;
    ILUp_precond(*sparILU, A->n, F, F);

    return 0;

}

/*----------------------------------------------------------------------------
 * ILUP preconditioner
 * Incomplete LU factorization with level of fill dropping
 *--------------------------------------------------------------------------*/
void ILUp(SparMAT* csmat, SparILU* lu, int p) {
    int n = csmat->n;
    int* jw, i, j, k, col, jpos, jrow;
    SparMAT* L;
    SparMAT* U;
    double* D;

    L = lu->L;
    U = lu->U;
    D = lu->D;

    jw = lu->work;
    /* set indicator array jw to -1 */
    for (j = 0; j < n; j++) jw[j] = -1;

    /* beginning of main loop */
    for (i = 0; i < n; i++) {
        /* set up the i-th row accroding to the nonzero
         * information from symbolic factorization */
        SPARILU_row(lu, i);

        /* setup array jw[], and initial i-th row */

        for (j = 0; j < L->nzcount[i]; j++) { /* initialize L part   */
            col = L->ja[i][j];
            jw[col] = j;
            L->ma[i][j] = 0;
        }

        jw[i] = i;
        D[i] = 0; /* initialize diagonal */

        for (j = 0; j < U->nzcount[i]; j++) { /* initialize U part   */
            col = U->ja[i][j];
            jw[col] = j;
            U->ma[i][j] = 0;
        }
        /* copy row from csmat into lu */
        for (j = 0; j < csmat->nzcount[i]; j++) {
            col = csmat->ja[i][j];
            jpos = jw[col];
            if (col < i)
                L->ma[i][jpos] = csmat->ma[i][j];
            else if (col == i)
                D[i] = csmat->ma[i][j];
            else
                U->ma[i][jpos] = csmat->ma[i][j];
        }
        /* eliminate previous rows */
        for (j = 0; j < L->nzcount[i]; j++) {
            jrow = L->ja[i][j];
            /* get the multiplier for row to be eliminated (jrow) */
            L->ma[i][j] *= D[jrow];

            /* combine current row and row jrow */
            for (k = 0; k < U->nzcount[jrow]; k++) {
                col = U->ja[jrow][k];
                jpos = jw[col];
                if (jpos == -1) continue;
                if (col < i)
                    L->ma[i][jpos] -= L->ma[i][j] * U->ma[jrow][k];
                else if (col == i)
                    D[i] -= L->ma[i][j] * U->ma[jrow][k];
                else
                    U->ma[i][jpos] -= L->ma[i][j] * U->ma[jrow][k];
            }
        }

        /* reset double-pointer to -1 ( U-part) */
        for (j = 0; j < L->nzcount[i]; j++) {
            col = L->ja[i][j];
            jw[col] = -1;
        }
        jw[i] = -1;
        for (j = 0; j < U->nzcount[i]; j++) {
            col = U->ja[i][j];
            jw[col] = -1;
        }

        if (D[i] == 0) {
            for (j = i + 1; j < n; j++) {
                L->ma[j] = NULL;
                U->ma[j] = NULL;
            }
            printf("fatal error: Zero diagonal found...\n");
            exit(1);
        }
        D[i] = 1.0 / D[i];
    }
    return;
}

