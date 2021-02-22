# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <assert.h>
# include "matrix.h"
# include "list.h"
# include "util.h"
# include "solvers.h"
# include "getTime.h"
# include "kcycle.h"
# include "dpa.h"
# include "error_correction.h"

void fill_diag (double *diag, int n, int ml) {
	int i;
	for (i=0; i<n; i++) diag[i] = 6.0;
	if (ml) { for (i=1; i<=n; i++) diag[i-1] += i*1.0e-2; }
};

matrix * generate_mat (int n, double *diag) {
	int i, j, n2=n*n, N=n2*n, ind_x, ind_y, ind_z;
	matrix *A = create_m(N, N, 7*N-6*n2);
	for (A->row_ptr[0]=0, j=0, i=0; i<N; i++) {
		ind_x = (i%n2)%n;
		ind_y = (i%n2)/n;
		ind_z = i/n2;
		if (ind_z>0) { A->col_ind[j] = i-n2; A->val[j] = -1.0; j++; }
		if (ind_y>0) { A->col_ind[j] =  i-n; A->val[j] = -1.0; j++; }
		if (ind_x>0) { A->col_ind[j] =  i-1; A->val[j] = -1.0; j++; }
		A->col_ind[j] = i; A->val[j] = A->diag[i] = diag[ind_x]; j++;
		if (ind_x<n-1) { A->col_ind[j] =  i+1; A->val[j] = -1.0; j++; }
		if (ind_y<n-1) { A->col_ind[j] =  i+n; A->val[j] = -1.0; j++; }
		if (ind_z<n-1) { A->col_ind[j] = i+n2; A->val[j] = -1.0; j++; }
		A->row_ptr[i+1] = j;
	}
	return A;
};

void RHS_vector (double *f, int n, double *diag) {
	int i, n2=n*n, N=n2*n, ind_x, ind_y, ind_z;
	for (i=0; i<N; i++) {
		ind_x = (i%n2)%n;
		ind_y = (i%n2)/n;
		ind_z = i/n2;
		f[i] = diag[ind_x];
		if (ind_z>0) { f[i] -= 1.0; }
		if (ind_y>0) { f[i] -= 1.0; }
		if (ind_x>0) { f[i] -= 1.0; }
		if (ind_x<n-1) { f[i] -= 1.0; }
		if (ind_y<n-1) { f[i] -= 1.0; }
		if (ind_z<n-1) { f[i] -= 1.0; }
	}
};

void read_matrix (char *filename, matrix **A, double **f, int symm) {
	int i, N, nnz, row, col, symm_nnz=0, nz=0;
	double val;
	char trash[256];
	Node *aux;
	list **row_elem;
	FILE *file = fopen(filename, "r");
	while (fscanf(file, "%%%[^\n]\n", trash)==1) {}
	assert(fscanf(file, "%d %*d %d\n", &N, &nnz)==2);
	row_elem = (list **) malloc(N*sizeof(list *)); for (i=0; i<N; i++) row_elem[i] = create_l();
	for (i=0; i<nnz; i++) {
		assert(fscanf(file, "%d %d %lf\n", &col, &row, &val)==3); // reading transpose
		if (val!=0.0) {
			insert_l_sort(row_elem[row-1], col-1, val);
			if (symm && row!=col) { symm_nnz++; insert_l_sort(row_elem[col-1], row-1, val); }
		} else nz++;
	}
	fclose(file);
	nnz += symm_nnz-nz;
	*A = create_m(N, N, nnz);
	*f = (double *) calloc((*A)->m, sizeof(double));
	for ((*A)->row_ptr[0]=0, i=0, row=0; row<N; row++) {
		for (aux=row_elem[row]->head; aux; aux=aux->next, i++) {
			(*A)->val[i] = aux->val;
			(*A)->col_ind[i] = aux->elem; if (aux->elem==row) (*A)->diag[row] = aux->val;
			(*f)[row] += aux->val;
		}
		(*A)->row_ptr[row+1] = i;
	}
	for (i=0; i<N; i++) destroy_l(row_elem[i]); free(row_elem);
};

void print_matrix (char *filename, matrix *A) {
	int i, j;
	FILE *file = fopen(filename, "w");
	fprintf(file, "%d %d %d\n", A->n, A->m, A->nnz);
	for (i=0; i<A->m; i++) {
		for (j=A->row_ptr[i]; j<A->row_ptr[i+1]; j++) fprintf(file, "%d %d %.6f\n", A->col_ind[j]+1, i+1, A->val[j]); // writing transpose
	}
	fclose(file);
};

void read_F(double *F) {
    int n;
    FILE *fp = fopen("mtx/matrix.F", "r");
    assert(fscanf(fp,"%d",&n)==1);
    for(int i=0; i<n; i++)
	assert(fscanf(fp,"%lf",&F[i])==1);
    fclose(fp);
}

extern double solver;
extern double setup;
extern double ciclov;
extern double relax;
extern double matching;
extern double galerk;
extern double transf;

int main (int argc, char **argv) {
	
	if (argc<9) {
		printf("\n  Usage: ./program <file/n> <symm/ml> <method> <tol> <lmax> <k> <omega> <precond> <NCL> <str_thr> <aggress> <refinePass> <truncq> <trunc_fact> <nr> <k0> <kw>\n");
		printf("  - <file/n> is either matrix filename (.mtx) for general matrix problem\n");
		printf("                or number of grid points, in each direction, for stencil matrix problem\n");
		printf("  - <symm/ml> is: (1) symmetric matrix [general matrix problem] / multilayer case [stencil matrix problem];\n");
		printf("                  (0) otherwise\n");
		printf("  - <method> is: (1) SOR[<omega>,<tol>,<lmax>];\n");
		printf("                 (2) CG[<tol>,<lmax>];\n");
		printf("                 (3) GMRES[<k>,<tol>,<lmax>,<precond>,<NCL>,<str_thr>,<aggress>,<refinePass>,<truncq>,<trunc_fact>,<omega>,<nr>,<k0>,<kw>];\n");
		printf("                 (4) LCD[<k>,<tol>,<lmax>,<precond>,<NCL>,<str_thr>,<aggress>,<refinePass>,<truncq>,<trunc_fact>,<omega>,<nr>];\n");
		printf("                 (5) AMG[<tol>,<lmax>,<NCL>,<str_thr>,<aggress>,<refinePass>,<truncq>,<trunc_fact>,<omega>,<nr>]\n");
		printf("                 (6) AMG K-Cycle[<k>,<tol>,<lmax>,<NCL>,<str_thr>,<aggress>,<refinePass>,<truncq>,<trunc_fact>,<omega>,<nr>,<k0>,<kw>]\n");
		printf("                 (7) DPA AMG[<tol>,<lmax>,<omega>,<NCL>,<str_thr>,<nr>]\n");
		printf("                 (8) Error Correction [<tol>,<lmax>,<omega>,<NCL>,<nr>]\n");
		printf("  - <tol> is desired tolerance\n");
		printf("  - <lmax> is maximum number of iterations\n");
		printf("  - <k> is number of vectors in Krylov basis\n");
		printf("  - <omega> is real parameter [in (0,2)] for SOR relaxation\n");
		printf("  - <precond> is: (0) none;\n");
		printf("                  (1) Gauss-Seidel;\n");
		printf("                  (2) AMG\n");
		printf("                  (3) DPA AMG\n");
		printf("                  (4) ILU p\n");
		printf("                  (5) AMG K-Ciclo\n");
		printf("                  (6) DPA K-Ciclo\n");
		printf("  - <NCL> is number of coarse-grid levels in AMG\n");
		printf("  - <str_thr> is real parameter [in (0,1)] for determining strong connections during setup\n");
		printf("                  for PREIS matching, use str_thr = -1\n");
		printf("  - <aggress> is 1 to apply aggressive coarsening on top level, 0 otherwise\n");
		printf("  - <refinePass> is 1 to apply refinement coarsening pass, 0 otherwise\n");
		printf("  - <truncq> is 1 to apply interpolation truncation, 0 otherwise\n");
		printf("  - <trunc_fact> is real parameter [in (0,1)] for truncation\n");
		printf("  - <nr> is number of relaxations performed per level (going down/up) during V-cycle solve\n");
		printf("  Example: ./program test.mtx 0 3 1e-8 1000 20 1.5 2 2 0.25 0 0 0 0 2\n\n");
		return 0;
	}
	int n, converg, symm=atoi(argv[2]), method=atoi(argv[3]), k, lmax, precond, NCL, aggress, refinePass, truncq, nr, k0, kw, ILU_p;
	double *diag=NULL, *f, *u, tol, omega, str_thr, trunc_fact, beta;
	// parametros dos metodos do marcelo
	if(method < 6) {
		tol		= atof(argv[4]);
		k		= atoi(argv[6]);
		lmax	= atoi(argv[5]);
		omega	= atof(argv[7]);
		precond	= atoi(argv[8]);
		NCL		= atoi(argv[9]);
		if (precond == 4) // ILUp
			ILU_p       = atoi(argv[9]);
		else {
			str_thr		= atof(argv[10]);
			aggress		= atoi(argv[11]);
			refinePass	= atoi(argv[12]);
			truncq		= atoi(argv[13]);
			trunc_fact	= atof(argv[14]);
			nr			= atoi(argv[15]);
		}
		if ((precond == 5) || (precond == 6)) {
			k0		= atof(argv[16]);
			kw		= atoi(argv[17]);
		}
	} else
	// parametros do k-cycle
	if (method == 6) {
		tol			= atof(argv[4]);
		k			= atoi(argv[6]);
		lmax		= atoi(argv[5]);
		omega		= atof(argv[7]);
		precond		= atoi(argv[8]);
		NCL			= atoi(argv[9]);
		str_thr		= atof(argv[10]);
		aggress		= atoi(argv[11]);
		refinePass	= atoi(argv[12]);
		truncq		= atoi(argv[13]);
		trunc_fact	= atof(argv[14]);
		nr			= atoi(argv[15]);
		k0			= atoi(argv[16]);
		kw			= atoi(argv[17]);
	} else
	if ((method == 7) || (method == 8)) {
		tol			= atof(argv[4]);
		lmax		= atoi(argv[5]);
		omega		= atof(argv[6]);
		NCL			= atoi(argv[7]);
		beta		= atof(argv[8]);
		nr			= atoi(argv[9]);
	}
	
	char *arg1;
	matrix *A, *A_s;
	n = strtol(argv[1], &arg1, 10);
	if (strcmp(arg1, argv[1])) { // stencil matrix
		printf("stencil_%d", n);
		diag = (double *) malloc(n*sizeof(double)); fill_diag(diag, n, symm);
		A = generate_mat(n, diag); //print_matrix("exM.mtx", A);
		f = (double *) malloc((n*n*n)*sizeof(double)); RHS_vector(f, n, diag);
	} else {
		printf("%s\n", arg1);
		read_matrix(arg1, &A, &f, symm);
	}
        //read_F(f);

	u = (double *) calloc(A->n, sizeof(double));
	switch (method) {
            case 1:  converg = SOR(A, f, u, omega, tol, lmax); break;
            case 2:  converg = CG(A, f, u, tol, lmax); break;
            case 3:  converg = GMRES(A, f, u, k, tol, lmax, precond, NCL, str_thr, aggress, refinePass, truncq, trunc_fact, omega, nr, ILU_p, k0, kw); break;
            case 4:  converg = LCD(A, f, u, k, tol, lmax, precond, NCL, str_thr, aggress, refinePass, truncq, trunc_fact, omega, nr, ILU_p); break;
            case 6:  converg = AMG_kcycle(A, f, u, k, NCL, str_thr, aggress, refinePass, truncq, trunc_fact, omega, nr, tol, lmax, k0, kw);  break;
            case 7:  converg = dpa_AMG (A, f, u, NCL, omega, nr, tol, lmax, beta);  break;
            case 8:  converg = error(A, f, u, NCL, omega, nr, tol, lmax);  break;
            default: converg = AMG(A, f, u, NCL, str_thr, aggress, refinePass, truncq, trunc_fact, omega, nr, tol, lmax);
	}
	printf("\n\tmin: %.6f \tmax: %.6f \tconverg: %d\n\n", dmin(u, A->n), dmax(u, A->n), converg);
	/*printf("\n\nRESULTADO:\n");
	for(int i=0; i<100; i++)
		printf("v[%d] = %lf \t",i,u[i]);
	printf("\n\n");*/
	free(f); free(u);
	if (strcmp(arg1, argv[1])) free(diag);
	destroy_m(A);
	return 0;
};

