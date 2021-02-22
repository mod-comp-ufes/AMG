# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include "amg.h"
# include "list.h"
# include "fib_heap.h"
# include "util.h"
# include "solvers.h"
# include "getTime.h"

matrix * new_standard_coarsening (matrix *A, int *CF_split, double *max_aik, int refinePass) {
	int i, j, k, l, lambda, max_i, connect;
	matrix *A_t = transpose(A);
	Node1 *aux, *prev;
	list1 *fine_points=NULL;
	node_fh *max, **fib_n = (node_fh **) malloc(A->m*sizeof(node_fh *));
	fib_heap *fib = make_heap();
	for (i=0; i<A->m; i++) {
		lambda = 0;
		for (j=A_t->row_ptr[i]; j<A_t->row_ptr[i+1]; j++) {
			if (A_t->col_ind[j]!=i && fabs(A_t->val[j])>=max_aik[A_t->col_ind[j]]) {
				lambda++;
			}
		}
		fib_n[i] = insert_heap(fib, i, -lambda);
	}
	if (refinePass) fine_points = create_l1();
	while (fib->n) {
		max = extract_min_heap(fib);
		lambda = max->key; max_i = max->v;
		free(max);
		if (!lambda) {
			CF_split[max_i] = 0; if (refinePass) insert_l1(fine_points, max_i);
			break;
		}
		CF_split[max_i] = 1; // variable indexed by max_i becomes C-variable
		for (j=A_t->row_ptr[max_i]; j<A_t->row_ptr[max_i+1]; j++) {
			if (/*A_t->col_ind[j]!=max_i &&*/ fabs(A_t->val[j])>=max_aik[A_t->col_ind[j]] && CF_split[A_t->col_ind[j]]<0) {
				CF_split[A_t->col_ind[j]] = 0; // variable indexed by A_t->col_ind[j] becomes F-variable
				if (refinePass) insert_l1(fine_points, A_t->col_ind[j]);
				delete_heap(fib, fib_n[A_t->col_ind[j]]);
				for (i=A->row_ptr[A_t->col_ind[j]]; i<A->row_ptr[A_t->col_ind[j]+1]; i++) { // updating lambda w.r.t. new F-variable
					if (/*A->col_ind[i]!=A_t->col_ind[j] &&*/ fabs(A->val[i])>=max_aik[A_t->col_ind[j]] && CF_split[A->col_ind[i]]<0) {
						decrease_key_heap(fib, fib_n[A->col_ind[i]], fib_n[A->col_ind[i]]->key-1);
					}
				}
			}
		}
	}
	free(fib_n);
	while (fib->n) {
		max = extract_min_heap(fib);
		CF_split[max->v] = 0; if (refinePass) insert_l1(fine_points, max->v);
		free(max);
	}
	free(fib);
	if (refinePass) {
		for (prev=NULL, aux=fine_points->head; aux;) {
			i = aux->elem;
			for (j=A_t->row_ptr[i]; j<A_t->row_ptr[i+1]; j++) {
				if (A_t->col_ind[j]!=i && !CF_split[A_t->col_ind[j]] && fabs(A_t->val[j])>=max_aik[A_t->col_ind[j]]) {
					connect = 0;
					for (k=A_t->row_ptr[A_t->col_ind[j]]; k<A_t->row_ptr[A_t->col_ind[j]+1]; k++) {
						if (CF_split[A_t->col_ind[k]] && fabs(A_t->val[k])>=max_aik[A_t->col_ind[k]]) {
							for (l=A->row_ptr[A_t->col_ind[k]]; l<A->row_ptr[A_t->col_ind[k]+1]; l++) {
								if (A->col_ind[l]==i && fabs(A->val[l])>=max_aik[A_t->col_ind[k]]) {
									connect = 1;
									break;
								}
							}
							if (connect) break;
						}
					}
					if (!connect) {
						CF_split[i] = 1;
						break;
					}
				}
			}
			if (CF_split[i]) { // update fine_points
				if (!prev) {
					fine_points->head = aux->next;
					free(aux); fine_points->n--;
					aux = fine_points->head;
				} else {
					prev->next = aux->next;
					free(aux); fine_points->n--;
					aux = prev->next;
				}
			} else { prev = aux; aux = aux->next; }
		}
		destroy_l1(fine_points);
	}
	return A_t;
};

void new_agressive_coarsening (matrix *A, int *CF_split, double *max_aik) {
	int i, j, k, l, lambda, max_i, *marked;
	matrix *A_t = new_standard_coarsening(A, CF_split, max_aik, 0);
	marked = (int *) calloc(A->m, sizeof(int));
	node_fh *max, **fib_n = (node_fh **) malloc(A->m*sizeof(node_fh *));
	fib_heap *fib = make_heap();
	for (i=0; i<A->m; i++) {
		if (CF_split[i]) {
			CF_split[i] = -1;
			lambda = 0;
			for (j=A_t->row_ptr[i]; j<A_t->row_ptr[i+1]; j++) {
				if (A_t->col_ind[j]!=i && fabs(A_t->val[j])>=max_aik[A_t->col_ind[j]]) {
					if (CF_split[A_t->col_ind[j]] && marked[A_t->col_ind[j]]!=i+1) {
						lambda++;
						marked[A_t->col_ind[j]] = i+1;
					}
					for (k=A_t->row_ptr[A_t->col_ind[j]]; k<A_t->row_ptr[A_t->col_ind[j]+1]; k++) {
						if (A_t->col_ind[k]!=A_t->col_ind[j] && CF_split[A_t->col_ind[k]] && fabs(A_t->val[k])>=max_aik[A_t->col_ind[k]] && marked[A_t->col_ind[k]]!=i+1) {
							lambda++;
							marked[A_t->col_ind[k]] = i+1;
						}
					}
				}
			}
			fib_n[i] = insert_heap(fib, i, -lambda);
		}
	}
	while (fib->n) {
		max = extract_min_heap(fib);
		lambda = max->key; max_i = max->v;
		free(max);
		if (!lambda) {
			CF_split[max_i] = 1;
			break;
		}
		CF_split[max_i] = 1; // variable indexed by max_i becomes C-variable
		for (j=A_t->row_ptr[max_i]; j<A_t->row_ptr[max_i+1]; j++) {
			if (A_t->col_ind[j]!=max_i && fabs(A_t->val[j])>=max_aik[A_t->col_ind[j]]) {
				if (CF_split[A_t->col_ind[j]]<0) {
					CF_split[A_t->col_ind[j]] = 0; // variable indexed by A_t->col_ind[j] becomes F-variable
					delete_heap(fib, fib_n[A_t->col_ind[j]]);
					for (i=A->row_ptr[A_t->col_ind[j]]; i<A->row_ptr[A_t->col_ind[j]+1]; i++) { // updating lambda w.r.t. new F-variable
						if (A->col_ind[i]!=A_t->col_ind[j] && fabs(A->val[i])>=max_aik[A_t->col_ind[j]]) {
							if (CF_split[A->col_ind[i]]<0 && marked[A->col_ind[i]]!=-(A_t->col_ind[j]+1)) {
								decrease_key_heap(fib, fib_n[A->col_ind[i]], fib_n[A->col_ind[i]]->key-1);
								marked[A->col_ind[i]] = -(A_t->col_ind[j]+1);
							}
							for (k=A->row_ptr[A->col_ind[i]]; k<A->row_ptr[A->col_ind[i]+1]; k++) {
								if (/*A->col_ind[k]!=A->col_ind[i] &&*/ fabs(A->val[k])>=max_aik[A->col_ind[i]] && CF_split[A->col_ind[k]]<0 && marked[A->col_ind[k]]!=-(A_t->col_ind[j]+1)) {
									decrease_key_heap(fib, fib_n[A->col_ind[k]], fib_n[A->col_ind[k]]->key-1);
									marked[A->col_ind[k]] = -(A_t->col_ind[j]+1);
								}
							}
						}
					}
				}
				for (k=A_t->row_ptr[A_t->col_ind[j]]; k<A_t->row_ptr[A_t->col_ind[j]+1]; k++) {
					if (/*A_t->col_ind[k]!=A_t->col_ind[j] &&*/ fabs(A_t->val[k])>=max_aik[A_t->col_ind[k]] && CF_split[A_t->col_ind[k]]<0) {
						CF_split[A_t->col_ind[k]] = 0; // variable indexed by A_t->col_ind[k] becomes F-variable
						delete_heap(fib, fib_n[A_t->col_ind[k]]);
						for (i=A->row_ptr[A_t->col_ind[k]]; i<A->row_ptr[A_t->col_ind[k]+1]; i++) { // updating lambda w.r.t. new F-variable
							if (A->col_ind[i]!=A_t->col_ind[k] && fabs(A->val[i])>=max_aik[A_t->col_ind[k]]) {
								if (CF_split[A->col_ind[i]]<0 && marked[A->col_ind[i]]!=-(A_t->col_ind[k]+1)) {
									decrease_key_heap(fib, fib_n[A->col_ind[i]], fib_n[A->col_ind[i]]->key-1);
									marked[A->col_ind[i]] = -(A_t->col_ind[k]+1);
								}
								for (l=A->row_ptr[A->col_ind[i]]; l<A->row_ptr[A->col_ind[i]+1]; l++) {
									if (/*A->col_ind[l]!=A->col_ind[i] &&*/ fabs(A->val[l])>=max_aik[A->col_ind[i]] && CF_split[A->col_ind[l]]<0 && marked[A->col_ind[l]]!=-(A_t->col_ind[k]+1)) {
										decrease_key_heap(fib, fib_n[A->col_ind[l]], fib_n[A->col_ind[l]]->key-1);
										marked[A->col_ind[l]] = -(A_t->col_ind[k]+1);
									}
								}
							}
						}
					}
				}
			}
		}
	}
	free(fib_n);
	while (fib->n) { max = extract_min_heap(fib); CF_split[max->v] = 1; free(max); } free(fib);
	free(marked);
	destroy_m(A_t);
};

/*int truncation_interp (list *row_elem, double trunc_fact) {
	int i=0;
	double max, sum1, sum2, rescale;
	Node *aux, *prev;
	max = sum1 = 0.0;
	for (aux=row_elem->head; aux; aux=aux->next) { sum1 += aux->val; if (fabs(aux->val)>max) max = fabs(aux->val); }
	sum2 = 0.0;
	for (prev=NULL, aux=row_elem->head; aux;) {
		if (fabs(aux->val)<trunc_fact*max) {
			i++;
			if (!prev) {
				row_elem->head = aux->next;
				free(aux); row_elem->n--;
				aux = row_elem->head;
			} else {
				prev->next = aux->next;
				free(aux); row_elem->n--;
				aux = prev->next;
			}
		} else {
			sum2 += aux->val;
			prev = aux; aux = aux->next;
		}
	}
	rescale = sum1/sum2;
	for (aux=row_elem->head; aux; aux=aux->next) { aux->val *= rescale; }
	return i;
};*/

int sign (double x) {
	return (x>0.0)? 1 : ((x<0.0)? -1 : 0);
};

void extended_interpolation (matrix *A, int *CF_split, int n_c, matrix **I_cf, double *max_aik, int truncq, double trunc_fact) {
	int i, j, k, nnz, sign_d_k, *interpSet, exist, *xb;
	double d_i, sum, val, *x;
	Node *aux;
	list **row_elem = (list **) malloc(A->m*sizeof(list *));
	interpSet = (int *) calloc(A->m, sizeof(int));
	xb = (int *) calloc(n_c, sizeof(int));
	x = (double *) malloc(n_c*sizeof(double));
	for (nnz=0, i=0; i<A->m; i++) {
		row_elem[i] = create_l();
		if (CF_split[i]) {
			// C-variables
			insert_l_tail(row_elem[i], CF_split[i]-1, 1.0);
			nnz++;
		} else {
			// F-variables
			for (j=A->row_ptr[i]; j<A->row_ptr[i+1]; j++) {
				// loop for defining interpolatory variables
				if (A->col_ind[j]!=i && fabs(A->val[j])>=max_aik[i]) {
					if (CF_split[A->col_ind[j]]) {
						// A->col_ind[j] E C_i^s
						interpSet[A->col_ind[j]] = i+1;
					} else {
						// A->col_ind[j] E F_i^s
						for (k=A->row_ptr[A->col_ind[j]]; k<A->row_ptr[A->col_ind[j]+1]; k++) {
							if (CF_split[A->col_ind[k]] && fabs(A->val[k])>=max_aik[A->col_ind[j]]) {
								// A->col_ind[k] E C_(A->col_ind[j])^s
								interpSet[A->col_ind[k]] = i+1;
							}
						}
					}
				}
			}
			d_i = A->diag[i];
			for (j=A->row_ptr[i]; j<A->row_ptr[i+1]; j++) {
				if (A->col_ind[j]!=i) {
					if (fabs(A->val[j])<max_aik[i]) {
						// A->col_ind[j] E N_i^w
						if (interpSet[A->col_ind[j]]!=i+1) {
							// A->col_ind[j] E N_i^w \ C_i
							d_i += A->val[j];
						}
					} else {
						if (!CF_split[A->col_ind[j]]) {
							// A->col_ind[j] E F_i^s
							sign_d_k = sign(A->diag[A->col_ind[j]]);
							for (exist=0, sum=0.0, k=A->row_ptr[A->col_ind[j]]; k<A->row_ptr[A->col_ind[j]+1]; k++) {
								if (interpSet[A->col_ind[k]]==i+1 && sign(A->val[k])!=sign_d_k) {
									sum += A->val[k];
								} else if (A->col_ind[k]==i && sign(A->val[k])!=sign_d_k) {
									exist = 1;
									val = A->val[k];
									sum += val;
								}
							}
							if (exist) {
								d_i += (A->val[j]*val)/sum;
							}
						}
					}
				}
			}
			for (j=A->row_ptr[i]; j<A->row_ptr[i+1]; j++) {
				if (interpSet[A->col_ind[j]]==i+1) {
					insert_l_tail(row_elem[i], CF_split[A->col_ind[j]]-1, 0.0);
					nnz++;
					xb[CF_split[A->col_ind[j]]-1] = i+1;
					x[CF_split[A->col_ind[j]]-1] = A->val[j];
				}
			}
			for (k=A->row_ptr[i]; k<A->row_ptr[i+1]; k++) {
				if (A->col_ind[k]!=i && fabs(A->val[k])>=max_aik[i] && !CF_split[A->col_ind[k]]) {
					// A->col_ind[k] E F_i^s
					sign_d_k = sign(A->diag[A->col_ind[k]]);
					for (sum=0.0, j=A->row_ptr[A->col_ind[k]]; j<A->row_ptr[A->col_ind[k]+1]; j++) {
						if ((interpSet[A->col_ind[j]]==i+1 || A->col_ind[j]==i) && sign(A->val[j])!=sign_d_k) {
							sum += A->val[j];
						}
					}
					for (j=A->row_ptr[A->col_ind[k]]; j<A->row_ptr[A->col_ind[k]+1]; j++) {
						if (interpSet[A->col_ind[j]]==i+1 && sign(A->val[j])!=sign_d_k) {
							if (xb[CF_split[A->col_ind[j]]-1]!=i+1) {
								insert_l_tail(row_elem[i], CF_split[A->col_ind[j]]-1, 0.0);
								nnz++;
								xb[CF_split[A->col_ind[j]]-1] = i+1;
								x[CF_split[A->col_ind[j]]-1] = (A->val[k]*A->val[j])/sum;
							} else {
								x[CF_split[A->col_ind[j]]-1] += (A->val[k]*A->val[j])/sum;
							}
						}
					}
				}
			}
			for (aux=row_elem[i]->head; aux; aux=aux->next) {
				aux->val = -x[aux->elem]/d_i;
			}
			//if (truncq) nnz -= truncation_interp(row_elem[i], trunc_fact); // 0.2
		}
	}
	free(interpSet);
	free(xb);
	free(x);
	*I_cf = create_m(A->m, n_c, nnz);
	for ((*I_cf)->row_ptr[0]=0, j=0, i=0; i<A->m; i++) {
		for (aux=row_elem[i]->head; aux; aux=aux->next, j++) {
			(*I_cf)->val[j] = aux->val;
			(*I_cf)->col_ind[j] = aux->elem;
			if (aux->elem==i) {
				(*I_cf)->diag[i] = aux->val;
			}
		}
		(*I_cf)->row_ptr[i+1] = j;
	}
	for (i=0; i<A->m; i++) {
		destroy_l(row_elem[i]);
	}
	free(row_elem);
};

void setupCL (matrix *A, int *CF_split, matrix **I_cf, matrix **I_cf_t, double str_thr, int aggress, int refinePass, int truncq, double trunc_fact, int level) {
	int i, j, n_c;
	double *max_aik;
	matrix *A_t;
	max_aik = (double *) calloc(A->m, sizeof(double));
	for (i=0; i<A->m; i++) {
		for (j=A->row_ptr[i]; j<A->row_ptr[i+1]; j++) {
			if (A->col_ind[j]!=i && fabs(A->val[j])>max_aik[i]) max_aik[i] = fabs(A->val[j]);
		}
		max_aik[i] *= str_thr;
	}
	if (level || !aggress) {
		A_t = new_standard_coarsening(A, CF_split, max_aik, refinePass);
		destroy_m(A_t);
	} else new_agressive_coarsening(A, CF_split, max_aik);
	for (n_c=0, i=0; i<A->m; i++) {
		if (CF_split[i]) { // i E C
			n_c++;
			CF_split[i] = n_c;
		}
	}
	extended_interpolation(A, CF_split, n_c, I_cf, max_aik, truncq, trunc_fact);
	free(max_aik);
	*I_cf_t = transpose(*I_cf);
};

matrix * coarse_matrix (matrix *A, matrix *I_cf, matrix *I_cf_t) { // A_c = (I_cf^t)A(I_cf)
	matrix *B, *C;
	B = matmat_sparseproduct(I_cf_t, A); 
	C = matmat_sparseproduct(B, I_cf);
	destroy_m(B);
	return C;
};

double AMG_setup (matrix *A_o, double *f_o, double *u_o, int NCL, double str_thr, int aggress, int refinePass, int truncq, double trunc_fact, matrix ***A, double ***f, double ***u, double **r, matrix ***I_cf, matrix ***I_cf_t) {
	int i, j, *CF_split, fix=1;
	double t;
	if (NCL==0) { fix = 0; NCL = 10; }
	*A       = (matrix **) malloc((NCL+1)*sizeof(matrix *));
	*f       = (double **) malloc((NCL+1)*sizeof(double *));
	CF_split = (int *)     malloc(A_o->m*sizeof(int));
	*I_cf    = (matrix **) malloc(NCL*sizeof(matrix *));
	*I_cf_t  = (matrix **) malloc(NCL*sizeof(matrix *));
	*u       = (double **) malloc((NCL+1)*sizeof(double *));
	(*A)[0] = A_o; (*f)[0] = f_o; (*u)[0] = u_o;
	// SETUP PHASE
	t = get_time();
	if (!fix) {
		NCL = 0;
		for (i=0; (*A)[i]->m>10000 && i<10; i++) {
			NCL++;
			for (j=0; j<(*A)[i]->m; j++) CF_split[j] = -1;
			setupCL((*A)[i], CF_split, &((*I_cf)[i]), &((*I_cf_t)[i]), str_thr, aggress, refinePass, truncq, trunc_fact, i);
                        (*A)[i+1] = coarse_matrix((*A)[i], (*I_cf)[i], (*I_cf_t)[i]);
		}
	} else {
		for (i=0; i<NCL; i++) {
			for (j=0; j<(*A)[i]->m; j++) CF_split[j] = -1;
			setupCL((*A)[i], CF_split, &((*I_cf)[i]), &((*I_cf_t)[i]), str_thr, aggress, refinePass, truncq, trunc_fact, i);
                        (*A)[i+1] = coarse_matrix((*A)[i], (*I_cf)[i], (*I_cf_t)[i]);
		}
	}
	t = (get_time() - t)/100.0;
	free(CF_split);
	for (i=1; i<NCL+1; i++) {
		(*f)[i] = (double *) malloc((*A)[i]->m*sizeof(double));
		(*u)[i] = (double *) malloc((*A)[i]->n*sizeof(double));
	}
	*r = (double *) malloc(A_o->m*sizeof(double));
	return t;
};

void AMG_destroy (int NCL, matrix **A, double **f, double **u, double *r, matrix **I_cf, matrix **I_cf_t) {
	int i;
	for (i=1; i<NCL+1; i++) { free(f[i]); free(u[i]); }
	free(f); free(u); free(r);
	for (i=1; i<NCL+1; i++) destroy_m(A[i]);
	for (i=0; i<NCL; i++) { destroy_m(I_cf[i]); destroy_m(I_cf_t[i]); }
	free(A); free(I_cf); free(I_cf_t);
};

void AMG_Vcycle (int NCL, matrix **A, double **f, double **u, double *r, matrix **I_cf, matrix **I_cf_t, double omega, int nr) {

	int i, j;
	for (i=0; i<NCL; i++) {
		SOR_relax(A[i], f[i], u[i], omega, nr);
		residual(r, f[i], A[i], u[i], 0, 0, 0);
		mat_vec(f[i+1], I_cf_t[i], r);
		for (j=0; j<A[i+1]->n; j++) u[i+1][j] = 0.0;
	}
	SOR_relax(A[i], f[i], u[i], omega, 2*nr);
	for (i--; i>=0; i--) {
		mat_vec_plus(u[i], I_cf[i], u[i+1]);
		SOR_relax(A[i], f[i], u[i], omega, nr);
	}
};

extern double transf;

void AMG_Vcycle_precond (int NCL, matrix **A, double **f, double **u, double *r, matrix **I_cf, matrix **I_cf_t, double omega, int nr) {
	int i, j;
	for (i=0; i<A[0]->m; i++) u[0][i] = 0.0; // initial guess
	for (i=0; i<NCL; i++) {
		SOR_relax(A[i], f[i], u[i], omega, nr);
		residual(r, f[i], A[i], u[i], 0, 0, 0);
		mat_vec(f[i+1], I_cf_t[i], r);
		for (j=0; j<A[i+1]->n; j++) u[i+1][j] = 0.0;
	}
	SOR_relax(A[i], f[i], u[i], omega, nr);
	SOR_relax_rev(A[i], f[i], u[i], omega, nr);
	for (i--; i>=0; i--) {
		mat_vec_plus(u[i], I_cf[i], u[i+1]);
		SOR_relax_rev(A[i], f[i], u[i], omega, nr);
	}
};

int AMG (matrix *A_o, double *f_o, double *u_o, int NCL, double str_thr, int aggress, int refinePass, int truncq, double trunc_fact, double omega, int nr, double tol, int lmax) {
	int i, k;
	double **f, **u, *r, delta, t, norm_f;
	matrix **A, **I_cf, **I_cf_t;
	t = AMG_setup(A_o, f_o, u_o, NCL, str_thr, aggress, refinePass, truncq, trunc_fact, &A, &f, &u, &r, &I_cf, &I_cf_t);	
	printf("	AMG(%d,%.2f,%d,%d,%d,%.2f,%.2f,%d,%e)", NCL, str_thr, aggress, refinePass, truncq, trunc_fact, omega, nr, tol);
	for (i=0; i<NCL+1; i++) printf("	%d	%d", A[i]->m, A[i]->nnz);
	printf("	%.6f", t);
	residual(r, f_o, A_o, u_o, 0, 0, 0);
	norm_f = norm_inf(f_o, A_o->m);
	delta = norm_inf(r, A_o->m)/norm_f;
	k = 0;
	t = get_time();
	while (k<lmax && delta>tol) {
		AMG_Vcycle(NCL, A, f, u, r, I_cf, I_cf_t, omega, nr);
		residual(r, f_o, A_o, u_o, 0, 0, 0);
		delta = norm_inf(r, A_o->m)/norm_f;
		k++;
		//printf("  %.6f\n", delta);
	}
	t = (get_time() - t)/100.0;
	printf("	%d	%.6f", k, t);
	AMG_destroy(NCL, A, f, u, r, I_cf, I_cf_t);
	return (delta<=tol)? 1: 0;
};

