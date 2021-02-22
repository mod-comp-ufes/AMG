#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include "getTime.h"
#include "list.h"
#include "matrix.h"
#include "solvers.h"
#include "util.h"
#include "counting.h"
#include "bcm_linmatch.h"

int pairs; // contabiliza quantas arestas tem o matching

int * matching_c(matrix * A, double beta/*, int checkDD*/) {
    matrix * A_t = transpose(A);
    A = sum_abs(A, A_t);
    destroy_m(A_t);

    int neg = 0;

    int i, j, max_j, max_j_exist, nc = 0, n = A->n;
    int * row_ptr = A->row_ptr, * col_ind = A->col_ind;
    double * val = A->val;
    Node* aux;

    int * rmatch = (int*) malloc((n + 1) * sizeof (int));

    counting * c = create_c(n);

    /*double * max;
    if (checkDD) {
            max = (double*)calloc(n,sizeof(double));
            for(i=0, line=0; i < nnz; i++) {
                    if (i == row_ptr[line+1]) line++;
                    if (val[i] > max[line]) max[line] = val[i];
            }
    }*/

    list ** s = (list**) malloc(n * sizeof (list*));
    for (i = 0; i < n; i++) s[i] = create_l();

    /*TODO: checkDD*/

    /*	Encontra maior coeficiente de cada variavel [max]	*/
    double * max;
    if (neg) {
        max = (double*) calloc(n, sizeof (double));
        for (i = 0; i < n; i++) {
            for (j = row_ptr[i]; j < row_ptr[i + 1]; j++)
                if (val[j] < 0) {
                    if (fabs(val[j]) > max[i]) {
                        /*printf("\n %f > %f",fabs(val[j]),max[i]);*/
                        max[i] = fabs(val[j]);
                    }
                    /*else
                        printf("\n %f <= %f",fabs(val[j]),max[i]);*/
                }
                /*else
                    printf("\n %f >= 0",val[j]);*/
        }
    } else {
        max = (double*) calloc(n, sizeof (double));
        for (i = 0; i < n; i++) {
            for (j = row_ptr[i]; j < row_ptr[i + 1]; j++)
                if (fabs(val[j]) > max[i]) {
                    max[i] = fabs(val[j]);
                }
        }
    }

    /*	Insere vizinhos fortemente conectados em s	*/
    if (neg) {
        for (i = 0; i < n; i++) {
            for (j = row_ptr[i]; j < row_ptr[i + 1]; j++)
                if (i != col_ind[j])
                    if (val[j] < -beta * max[i]) {
                        //if (s[i]->n == 0) printf("\n [%d]",i);
                        insert_l_tail(s[i], 0, col_ind[j]);
                    }
        }
    } else {
        for (i = 0; i < n; i++) {
            for (j = row_ptr[i]; j < row_ptr[i + 1]; j++)
                if (i != col_ind[j])
                    if (val[j] > beta * max[i])
                        insert_l_tail(s[i], 0, col_ind[j]);
        }
    }

    /*	insere no na estrutura do counting sort	*/
    for (i = 0; i < n; i++) {
        insert_c(c, i, (s[i]->n)+(c->n));
    }
    int k = 0;
    while (c->min_m < (2 * n)) {
        k++;
        /* Passo 1 */
        i = c->count[c->min_m]->head->i;
        nc++;

        /* Passo 2 */
        max_j_exist = 0;
        max_j = row_ptr[i];
        for (j = row_ptr[i]; j < row_ptr[i + 1]; j++)
            if (i != col_ind[j])
                if (c->nodes[col_ind[j]]) {
                    if ((val[j] >= val[max_j]) || (max_j_exist==0)) {
                        max_j = j;
                    }
                    max_j_exist = 1;
                }
        j = col_ind[max_j];

        /* passo 3 */
        if (check_in_l(s[i], j) && max_j_exist) { // se j E Si
            /* Passo 4 */
            rmatch[i] = j;
            rmatch[j] = i;
            remove_c(c, i);
            remove_c(c, j);

            /* Passo 5 */
            for (aux = s[i]->head; aux; aux = aux->next) {
                if (c->nodes[(int) (aux->val)]) {
                    move_c(c, (int) (aux->val), c->m[(int) (aux->val)] - 1);
                }
            }
            for (aux = s[j]->head; aux; aux = aux->next) {
                if (c->nodes[(int) (aux->val)]) {
                    move_c(c, (int) (aux->val), c->m[(int) (aux->val)] - 1);
                }
            }
        } else {
            /* Passo 4 */
            rmatch[(int) (i)] = -1;
            remove_c(c, i);

            /* Passo 5 */
            for (aux = s[i]->head; aux; aux = aux->next) {
                if (c->nodes[(int) (aux->val)]) {
                    move_c(c, (int) (aux->val), c->m[(int) (aux->val)] - 1);
                }
            }
        }
    }
    
    free(max);
    for (i = 0; i < n; i++) destroy_l(s[i]);
    free(s);
    destroy_c(c);
    destroy_m(A);

    rmatch[n] = n - nc;

    return rmatch;
}

int * matching_c_beta0(matrix * A) {
    matrix * A_t = transpose(A);
    A = sum_abs(A, A_t);
    destroy_m(A_t);

    int i, j, max_j, max_j_exist, nc = 0, n = A->n;
    int * row_ptr = A->row_ptr, * col_ind = A->col_ind;
    double * val = A->val;

    int * rmatch = (int*) malloc((n + 1) * sizeof (int));

    counting * c = create_c(n);

    /*double * max;
    if (checkDD) {
            max = (double*)calloc(n,sizeof(double));
            for(i=0, line=0; i < nnz; i++) {
                    if (i == row_ptr[line+1]) line++;
                    if (val[i] > max[line]) max[line] = val[i];
            }
    }*/

    /*TODO: checkDD*/

    /*	insere no na estrutura do counting sort	*/
    for (i = 0; i < n; i++) {
        insert_c(c, i, (row_ptr[i+1] - row_ptr[i])+(c->n));
    }
    int k = 0;
    while (c->min_m < (2 * n)) {
        k++;
        /* Passo 1 */
        i = c->count[c->min_m]->head->i;
        nc++;

        /* Passo 2 */
        max_j_exist = 0;
        max_j = row_ptr[i];
        for (j = row_ptr[i]; j < row_ptr[i + 1]; j++)
            if (i != col_ind[j])
                if (c->nodes[col_ind[j]]) {
                    if ((val[j] >= val[max_j]) || (max_j_exist==0)) {
                        max_j = j;
                    }
                    max_j_exist = 1;
                }
        j = col_ind[max_j];

        /* passo 3 */
        if (max_j_exist) { // se j E Si
            /* Passo 4 */
            rmatch[i] = j;
            rmatch[j] = i;
            remove_c(c, i);
            remove_c(c, j);

            /* Passo 5 */
            for (int k = row_ptr[i]; k < row_ptr[i+1]; k++) {
                if (c->nodes[(int)(col_ind[k])]) {
                    move_c(c, (int)(col_ind[k]), c->m[(int) (col_ind[k])] - 1);
                }
            }
            for (int k = row_ptr[j]; k < row_ptr[j+1]; k++) {
                if (c->nodes[(int) (col_ind[k])]) {
                    move_c(c, (int)(col_ind[k]), c->m[(int) (col_ind[k])] - 1);
                }
            }
        } else {
            /* Passo 4 */
            rmatch[(int) (i)] = -1;
            remove_c(c, i);

            /* Passo 5 */
            for (int k = row_ptr[i]; k < row_ptr[i+1]; k++) {
                if (c->nodes[(int)(col_ind[k])]) {
                    move_c(c, (int)(col_ind[k]), c->m[(int) (col_ind[k])] - 1);
                }
            }
        }
    }
    
    destroy_m(A);

    rmatch[n] = n - nc;

    return rmatch;
}

matrix * dpa_galerkin(matrix * A, int * rmatch) {
    /*
            O operador de galerkin consiste no produto matriz*matriz*matriz P_t*A*P
		
            O prolongador P consiste numa matriz com altura igual a da matriz A e largura igual ao numero de agregados
            onde cada linha i tem exatamente um numero 1 na coluna j indicando que a variavel i pertence ao agregado j
		
            Uma vez que a matriz P tem exatamente um numero 1 por linha, e no maximo dois numeros 1 por coluna, nao ha
            necessidade de efetuar o produto matriz*matriz completo, que eh de ordem cubica.
		
            Esta funcao usa a seguinte estrategia: somam-se as linhas que fazem parte de um mesmo agregado, e em cada
            linha, somam-se os elementos que fazem parte de um mesmo agregado. O algoritmo soma elementos das linhas
            ao mesmo tempo em que soma as linhas, logo, a matriz eh lida apenas uma vez e sequencialmente, e portanto
            a funcao eh de ordem nnz (numero de nao nulos).
		
            Como nao e possivel prever o numero de elementos nao nulos da matriz resultante, primeiro a matriz eh
            criada numa lista de adjacencias, e depois transferida para CSR.
     */

    matrix * resulting_matrix;
    Node *aux;
    /* actual_line armazena qual linha da matriz resultante esta sendo preenchida atualmente */
    int actual_line = 0, i, nnz = 0, n = A->n, resulting_n = n - rmatch[n];
    /*
            line_1_ini e line_2_ini indicam os indices onde comecam as duas linhas a serem somadas, nos vetores val e
            col_ind do CSR da matriz original.
            Do mesmo modo, line_1_end e line_2_end indicam o termino e line_1_actual e line_2_actual indicam os
            elementos sendo lidos atualmente.		
     */
    int line_1_ini, line_1_end, line_1_actual, line_2_ini, line_2_end, line_2_actual;
    list **row_elem = (list **) malloc(resulting_n * sizeof (list *));
    for (i = 0; i < resulting_n; i++) row_elem[i] = create_l();
    /*
            Quando se insere um elemento tal que seu par nao foi inserido ainda, o ponteiro para o novo 'node'
            inserido eh guardado no vetor pair no indice do par. Quando o outro elemento do par for lido, em vez de
            inserir um novo 'node' o seu valor sera somado ao 'node' apontado no vetor.
     */
    Node **pair = (Node **) malloc(n * sizeof (Node *));
    /*
            Para saber se o par ja foi inserido ou nao, eh so verificar se pair_aux na posicao do 'node' contem o numero
            da linha. Apos inserir um elemento, armazena-se o numero da linha na posicao do 'node' no pair_aux para
            ajudar mais tarde.
     */
    int *pair_aux = (int*) malloc(n * sizeof (int));
    for (i = 0; i < n; i++) pair_aux[i] = -1;
    /*
            Pair_ind indica a qual agregado cada variavel pertence, o que eh necessario para saber o col_ind ao gerar
            o CSR da matriz resultante. 
     */
    int *pair_ind = (int*) malloc(n * sizeof (int));
    /*	atalhos		*/
    double * val = A->val;
    double * diag = A->diag;
    int * col_ind = A->col_ind;
    int * row_ptr = A->row_ptr;

    for (i = 0; i < n; i++) {
        /*
                Percorre as linhas da matriz original. Se a linha nao tiver par, eh transferida diretamente para
                a matriz resultante, pois nao ha nada a se somar a ela. Se ela tiver par, somam-se as duas linhas

                Como estamos lendo linha por linha, iremos passar pelas duas linhas de cada par, ou seja, cada par
                sera visitado duas vezes. Para processar cada par apenas uma vez, eh checado se o indice da linha
                eh menor que o indice do par. Se nao for, ignora. Assim cada par eh processado exatamente uma vez.
         */
        if (rmatch[i] == -1) {
            /*	rmatch[i] == -1 significa que i nao tem par	*/

            /*	calcula-se o inicio e o fim da linha a ser transferida		*/
            line_1_ini = line_1_actual = row_ptr[i];
            line_1_end = row_ptr[i + 1];

            while (line_1_actual < line_1_end) {
                /*	line_1_actual < line_1_end significa que a linha ainda nao acabou	*/
                if (pair_aux[col_ind[line_1_actual]] == i) {
                    /*	se o par ja foi inserido, soma o valor e limpa o vetor	*/
                    pair[col_ind[line_1_actual]]->val += val[line_1_actual];
                } else {
                    /*	se o par nao foi inserido ainda, insere novo 'node' e contabiliza	*/
                    /*if ((row_elem[actual_line]->tail) && (row_elem[actual_line]->tail->elem >= col_ind[line_1_actual])) {
                            printf("\na[%d][%d] = %f",i,col_ind[line_1_actual],val[line_1_actual]);
                            printf("\n");
                            for(int j=line_1_ini; j<line_1_end; j++) printf("%d ",col_ind[j]);
                            printf("\n");
                            exit(0);
                    }*/
                    if (rmatch[col_ind[line_1_actual]] == -1)
                        insert_l_tail(row_elem[actual_line], col_ind[line_1_actual], val[line_1_actual]);
                    else if (col_ind[line_1_actual] > rmatch[col_ind[line_1_actual]])
                        insert_l_sort(row_elem[actual_line], rmatch[col_ind[line_1_actual]], val[line_1_actual]);
                    else {
                        insert_l_tail(row_elem[actual_line], col_ind[line_1_actual], val[line_1_actual]);
                        /*	se houver par, grava o ponteiro para novo 'node' para posterior soma dos valores	*/
                        pair[rmatch[col_ind[line_1_actual]]] = row_elem[actual_line]->tail;
                        pair_aux[rmatch[col_ind[line_1_actual]]] = i;
                    }
                    nnz++;
                }
                /*	avanca na leitura da linha	*/
                line_1_actual++;
            }
            /*	o indice do elemento na nova matriz eh armazenada em pair_ind, para obter o col_ind mais tarde	*/
            pair_ind[i] = actual_line;
            actual_line++;
        } else
            if (i < rmatch[i]) {
            /*	i < rmatch[i] significa que o indice da linha i eh menor que o indice do par de i	*/

            /*	calcula-se inicio e fim das duas linhas a serem somadas		*/
            line_1_ini = line_1_actual = row_ptr[i];
            line_1_end = row_ptr[i + 1];
            line_2_ini = line_2_actual = row_ptr[rmatch[i]];
            line_2_end = row_ptr[rmatch[i] + 1];

            /*	enquanto houver coisa para transferir	*/
            while ((line_1_actual < line_1_end) || (line_2_actual < line_2_end)) {
                /*	Se a linha 2 acabou, certamente a linha 1 nao acabou, pois senao o while teria dado false.
                        Logo,  transfere a linha 1 sem somar nada pois nao ha nada a somar	*/
                if (line_2_actual == line_2_end) {
                    if (pair_aux[col_ind[line_1_actual]] == i) {
                        pair[col_ind[line_1_actual]]->val += val[line_1_actual];
                    } else {
                        /*if ((row_elem[actual_line]->tail) && (row_elem[actual_line]->tail->elem >= col_ind[line_1_actual])) {
                                printf("\nb[%d][%d] = %f",i,col_ind[line_1_actual],val[line_1_actual]);
                                printf("\n[%d] ",i);
                                for(int j=line_1_ini; j<line_1_end; j++) printf("%d ",col_ind[j]);
                                printf("\n[%d] ",rmatch[i]);
                                for(int j=line_2_ini; j<line_2_end; j++) printf("%d ",col_ind[j]);
                                printf("\n");
                                exit(0);
                        }*/
                        if (rmatch[col_ind[line_1_actual]] == -1)
                            insert_l_tail(row_elem[actual_line], col_ind[line_1_actual], val[line_1_actual]);
                        else if (col_ind[line_1_actual] > rmatch[col_ind[line_1_actual]])
                            insert_l_sort(row_elem[actual_line], rmatch[col_ind[line_1_actual]], val[line_1_actual]);
                        else {
                            insert_l_tail(row_elem[actual_line], col_ind[line_1_actual], val[line_1_actual]);
                            pair[rmatch[col_ind[line_1_actual]]] = row_elem[actual_line]->tail;
                            pair_aux[rmatch[col_ind[line_1_actual]]] = i;
                        }
                        nnz++;
                    }
                    line_1_actual++;
                } else
                    /*	Se a linha 2 nao acabou e a linha 1 sim, transfere a linha 2 sem somar nada	*/
                    /*	Se as duas linhas nao tiverem acabado, mas o indice dos elementos lidos forem diferentes,
                            transfere o de indice menor, a fim de manter a ordem	*/
                    if ((line_1_actual == line_1_end) || (col_ind[line_1_actual] > col_ind[line_2_actual])) {
                    if (pair_aux[col_ind[line_2_actual]] == i) {
                        pair[col_ind[line_2_actual]]->val += val[line_2_actual];
                    } else {
                        /*if ((row_elem[actual_line]->tail) && (row_elem[actual_line]->tail->elem >= col_ind[line_2_actual])) {
                                printf("\nc[%d][%d] = %f",i,col_ind[line_2_actual],val[line_2_actual]);
                                printf("\n[%d] ",i);
                                for(int j=line_1_ini; j<line_1_end; j++) printf("%d ",col_ind[j]);
                                printf("\n[%d] ",rmatch[i]);
                                for(int j=line_2_ini; j<line_2_end; j++) printf("%d ",col_ind[j]);
                                printf("\n");
                                exit(0);
                        }*/
                        if (rmatch[col_ind[line_2_actual]] == -1)
                            insert_l_tail(row_elem[actual_line], col_ind[line_2_actual], val[line_2_actual]);
                        else if (col_ind[line_2_actual] > rmatch[col_ind[line_2_actual]])
                            insert_l_sort(row_elem[actual_line], rmatch[col_ind[line_2_actual]], val[line_2_actual]);
                        else {
                            insert_l_tail(row_elem[actual_line], col_ind[line_2_actual], val[line_2_actual]);
                            pair[rmatch[col_ind[line_2_actual]]] = row_elem[actual_line]->tail;
                            pair_aux[rmatch[col_ind[line_2_actual]]] = i;
                        }
                        nnz++;
                    }
                    line_2_actual++;
                } else
                    if (col_ind[line_1_actual] < col_ind[line_2_actual]) {
                    if (pair_aux[col_ind[line_1_actual]] == i) {
                        pair[col_ind[line_1_actual]]->val += val[line_1_actual];
                    } else {
                        /*if ((row_elem[actual_line]->tail) && (row_elem[actual_line]->tail->elem >= col_ind[line_1_actual])) {
                                printf("\nd[%d][%d] = %f",i,col_ind[line_1_actual],val[line_1_actual]);
                                printf("\n[%d] ",i);
                                for(int j=line_1_ini; j<line_1_end; j++) printf("%d ",col_ind[j]);
                                printf("\n[%d] ",rmatch[i]);
                                for(int j=line_2_ini; j<line_2_end; j++) printf("%d ",col_ind[j]);
                                printf("\n");
                                exit(0);
                        }*/
                        if (rmatch[col_ind[line_1_actual]] == -1)
                            insert_l_tail(row_elem[actual_line], col_ind[line_1_actual], val[line_1_actual]);
                        else if (col_ind[line_1_actual] > rmatch[col_ind[line_1_actual]])
                            insert_l_sort(row_elem[actual_line], rmatch[col_ind[line_1_actual]], val[line_1_actual]);
                        else {
                            insert_l_tail(row_elem[actual_line], col_ind[line_1_actual], val[line_1_actual]);
                            pair[rmatch[col_ind[line_1_actual]]] = row_elem[actual_line]->tail;
                            pair_aux[rmatch[col_ind[line_1_actual]]] = i;
                        }
                        nnz++;
                    }
                    line_1_actual++;
                } else
                    /*	O unico caso que resta eh o que as duas linhas nao acabaram e os dois elementos lidos tem
                            o mesmo indice. Neste caso, transfere a soma. Assim, teremos somado as duas linhas	*/
                    //if (col_ind[line_1_actual] == col_ind[line_2_actual]) 
                {
                    if (pair_aux[col_ind[line_1_actual]] == i) {
                        pair[col_ind[line_1_actual]]->val += val[line_1_actual] + val[line_2_actual];
                    } else {
                        /*if ((row_elem[actual_line]->tail) && (row_elem[actual_line]->tail->elem >= col_ind[line_1_actual])) {
                                printf("\ne[%d][%d] = %f",i,col_ind[line_1_actual],val[line_1_actual]+val[line_2_actual]);
                                printf("\n[%d] ",i);
                                for(int j=line_1_ini; j<line_1_end; j++) printf("%d ",col_ind[j]);
                                printf("\n[%d] ",rmatch[i]);
                                for(int j=line_2_ini; j<line_2_end; j++) printf("%d ",col_ind[j]);
                                printf("\n");
                                exit(0);
                        }*/
                        if (rmatch[col_ind[line_1_actual]] == -1)
                            insert_l_tail(row_elem[actual_line], col_ind[line_1_actual], val[line_1_actual] + val[line_2_actual]);
                        else if (col_ind[line_1_actual] > rmatch[col_ind[line_1_actual]])
                            insert_l_sort(row_elem[actual_line], rmatch[col_ind[line_1_actual]], val[line_1_actual] + val[line_2_actual]);
                        else {
                            insert_l_tail(row_elem[actual_line], col_ind[line_1_actual], val[line_1_actual] + val[line_2_actual]);
                            pair[rmatch[col_ind[line_1_actual]]] = row_elem[actual_line]->tail;
                            pair_aux[rmatch[col_ind[line_1_actual]]] = i;
                        }
                        nnz++;
                    }
                    line_1_actual++;
                    line_2_actual++;
                }
            }
            /*	o indice do elemento na nova matriz eh armazenada em pair_ind, para obter o col_ind mais tarde	*/
            pair_ind[i] = pair_ind[rmatch[i]] = actual_line;
            actual_line++;
        }
    }

    resulting_matrix = create_m(resulting_n, resulting_n, nnz);
    val = resulting_matrix->val;
    diag = resulting_matrix->diag;
    col_ind = resulting_matrix->col_ind;
    row_ptr = resulting_matrix->row_ptr;
    int count = 0;
    row_ptr[0] = 0;
    for (i = 0; i < resulting_n; i++) {
        for (aux = row_elem[i]->head; aux; aux = aux->next, count++) {
            val[count] = aux->val;
            col_ind[count] = pair_ind[aux->elem];
            if (pair_ind[aux->elem] == i) diag[i] = aux->val;
            //if(i==1838)printf("\n[%d][%d/%d->%d]=%f",i,aux->elem,rmatch[aux->elem],pair_ind[aux->elem],aux->val);
            //printf("\n[%d][%d]=%f",i,pair_ind[aux->elem],aux->val); // descomente esta linha para printar a matriz
        }
        row_ptr[i + 1] = count;
    }

    for (i = 0; i < resulting_n; i++) destroy_l(row_elem[i]);
    free(row_elem);
    free(pair);
    free(pair_aux);
    free(pair_ind);

    return resulting_matrix;
}

void dpa_fine2coarse(int * match, double * fine, int n, double * coarse) {

    for (int i = 0, j = 0; i < n; i++)
        if (match[i] == -1)
            coarse[j++] = fine[i];
        else if (i < match[i])
            coarse[j++] = fine[i] + fine[match[i]];
}

void dpa_coarse2fine(int * match, double * coarse, int n, double * fine) {
    //double * fine = (double *)malloc(n*sizeof(double));
    for (int i = 0, j = 0; i < n; i++)
        if (match[i] == -1)
            fine[i] += coarse[j++];
        else if (i < match[i]) {
            fine[i] += coarse[j];
            fine[match[i]] += coarse[j++];
        }

    //return fine;
}

// EXPERIMENTAL. não use se não souber do que se trata
void dpa_coarse2fine_pure(int * match, double * coarse, int n, double * fine) {
    //double * fine = (double *)malloc(n*sizeof(double));
    for (int i = 0, j = 0; i < n; i++)
        if (match[i] == -1)
            fine[i] = coarse[j++];
        else if (i < match[i]) {
            fine[i] = coarse[j];
            fine[match[i]] = coarse[j++];
        }

    //return fine;
}

double dpa_setup(matrix *A_o, double *f_o, double *u_o, int NCL, matrix ***A, double ***f, double ***u, double **r, int ***match, double beta) {
    int i;
    double t, t2;
    NCL *= 2; // DPA
    *A = (matrix **) malloc((NCL + 1) * sizeof (matrix *));
    *f = (double **) malloc((NCL + 1) * sizeof (double *));
    *u = (double **) malloc((NCL + 1) * sizeof (double *));
    *match = (int **) malloc(NCL * sizeof (int *));
    (*A)[0] = A_o;
    (*f)[0] = f_o;
    (*u)[0] = u_o;
    matrix *A_s, *A_t;
    
    /*struct mc64_control mc_cntrl;
    struct mc64_info    mc_info;*/  
    
    // SETUP PHASE
    t = get_time();
    
    // gera A[i+1] com base em A[i]
    if (beta == -1) 
        for (i = 0; i < NCL; i++) {
            t2 = get_time();
            A_t = transpose((*A)[i]);
            A_s = sum_upper((*A)[i], A_t);
            (*match)[i] = bcm_CSRMatrixHMatch(A_s);
            (*A)[i + 1] = dpa_galerkin((*A)[i], (*match)[i]);
            //printf("\ni = %d \t n = %d \t m = %d \t nnz = %d \t rm = %d",i,(*A)[i]->n,(*A)[i]->m,(*A)[i]->nnz,(*match)[i][(*A)[i]->n]);
            destroy_m(A_t);
            destroy_m(A_s);
        }
    else
        for (i = 0; i < NCL; i++) {
            if (beta)
                (*match)[i] = matching_c((*A)[i], beta);
            else
                (*match)[i] = matching_c_beta0((*A)[i]);
            (*A)[i + 1] = dpa_galerkin((*A)[i], (*match)[i]);
        }

    for (i = 1; i < (NCL + 1); i++) {
        (*f)[i] = (double *) malloc((*A)[i]->m * sizeof (double));
        (*u)[i] = (double *) malloc((*A)[i]->n * sizeof (double));
    }
    *r = (double *) malloc(A_o->m * sizeof (double));
    
    t = (get_time() - t) / 100.0;
    
    return t;
};

void dpa_Vcycle_recursivo(int NCL, matrix **A, double **f, double **u, double *r, double omega, int nr, int i, int **match) {

    if (i < NCL) {

        SOR_relax(A[i], f[i], u[i], omega, nr);
        residual(r, f[i], A[i], u[i], 0, 0, NULL);
        //mat_vec(f[i+1], I_cf_t[i], r); // f[i+1] <- r
        dpa_fine2coarse(match[i], r, A[i]->n,f[i + 1]);

        for (int j = 0; j < A[i + 1]->n; j++) u[i + 1][j] = 0.0;

        dpa_Vcycle_recursivo(NCL, A, f, u, r, omega, nr, i + 1, match);

        //mat_vec_plus(u[i], I_cf[i], u[i+1]); // u[i] <- u[i+1]
        /*u[i] = */dpa_coarse2fine(match[i], u[i + 1], A[i]->n, u[i]);
        SOR_relax(A[i], f[i], u[i], omega, nr);
    } else {
        SOR_relax(A[i], f[i], u[i], omega, 2 * nr);
    }
}

void dpa_Vcycle(int NCL, matrix **A, double **f, double **u, double *r, double omega, int nr, int **match) {

    int i, j;
    for (i = 0; i < (2 * NCL); i++) {
        if ((i % 2) == 0)
            SOR_relax(A[i], f[i], u[i], omega, nr);
        residual(r, f[i], A[i], u[i], 0, 0, NULL);
        //mat_vec(f[i+1], I_cf_t[i], r); // f[i+1] <- r
        dpa_fine2coarse(match[i], r, A[i]->n,f[i + 1]);
        for (j = 0; j < A[i + 1]->n; j++) u[i + 1][j] = 0.0;
    }
    SOR_relax(A[i], f[i], u[i], omega, 2 * nr);
    for (i--; i >= 0; i--) {
        dpa_coarse2fine(match[i], u[i + 1], A[i]->n, u[i]);
        if ((i % 2) == 0)
            SOR_relax(A[i], f[i], u[i], omega, nr);
    }
};

// na descida, primeiro faz correcao por erro, depois NAO faz
void dpa_Vcycle2(int NCL, matrix **A, double **f, double **u, double *r, double omega, int nr, int **match) {

    int i, j;
    double t;
    NCL*=2;
    for (i = 0; i < NCL; i+=2) {
        SOR_relax(A[i], f[i], u[i], omega, nr);
        residual(r, f[i], A[i], u[i], 0, 0, NULL);
        dpa_fine2coarse(match[i],   r,      A[i]->n  , f[i+1]); // f[i+1] <- r
        dpa_fine2coarse(match[i+1], f[i+1], A[i+1]->n, f[i+2]); // f[i+2] <- f[i+1]
        for (j = 0; j < A[i+1]->n; j++) u[i+1][j] = 0.0;
        for (j = 0; j < A[i+2]->n; j++) u[i+2][j] = 0.0;
    }
    SOR_relax(A[i], f[i], u[i], omega, 2 * nr);
    for (; i > 0; i-=2) {
        dpa_coarse2fine_pure(match[i-1], u[i],   A[i-1]->n, u[i-1]); // u[i-1] <- u[i]
        dpa_coarse2fine     (match[i-2], u[i-1], A[i-2]->n, u[i-2]); // u[i-2] <- u[i-1]
        SOR_relax(A[i-2], f[i-2], u[i-2], omega, nr);
    }
};

// na descida, primeiro NAO faz correcao por erro, depois faz
void dpa_Vcycle3(int NCL, matrix **A, double **f, double **u, double *r, double omega, int nr, int **match) {

    int i, j;
    double t;
    NCL*=2;
    for (i = 0; i < NCL; i+=2) {
        SOR_relax(A[i], f[i], u[i], omega, nr);
        dpa_fine2coarse(match[i],   f[i], A[i]->n  , f[i+1]); // f[i+1] <- f[i]
        residual(r, f[i+1], A[i+1], u[i+1], 0, 0, NULL);
        dpa_fine2coarse(match[i+1], r,    A[i+1]->n, f[i+2]); // f[i+2] <- r
        for (j = 0; j < A[i+1]->n; j++) u[i+1][j] = 0.0;
        for (j = 0; j < A[i+2]->n; j++) u[i+2][j] = 0.0;
    }
    SOR_relax(A[i], f[i], u[i], omega, 2 * nr);
    for (; i > 0; i-=2) {
        dpa_coarse2fine     (match[i-1], u[i],   A[i-1]->n, u[i-1]); // u[i-1] <- u[i]
        dpa_coarse2fine_pure(match[i-2], u[i-1], A[i-2]->n, u[i-2]); // u[i-2] <- u[i-1]
        SOR_relax(A[i-2], f[i-2], u[i-2], omega, nr);
    }
};

void dpa_Vcycle_precond(int NCL, matrix **A, double **f, double **u, double *r, double omega, int nr, int **match) {

    int i, j;
    double t;
    for (i=0; i<A[0]->m; i++) u[0][i] = 0.0;
    for (i = 0; i < (2 * NCL); i++) {
        if ((i % 2) == 0)
            SOR_relax(A[i], f[i], u[i], omega, nr);
        residual(r, f[i], A[i], u[i], 0, 0, NULL);
        dpa_fine2coarse(match[i], r, A[i]->n,f[i + 1]);
        for (j = 0; j < A[i + 1]->n; j++) u[i + 1][j] = 0.0;
    }
    SOR_relax(A[i], f[i], u[i], omega, nr);
    SOR_relax_rev(A[i], f[i], u[i], omega, nr);
    for (i--; i >= 0; i--) {
        dpa_coarse2fine(match[i], u[i + 1], A[i]->n, u[i]);
        if ((i % 2) == 0)
            SOR_relax_rev(A[i], f[i], u[i], omega, nr);
    }
};

void dpa_Vcycle_precond2(int NCL, matrix **A, double **f, double **u, double *r, double omega, int nr, int **match) {

    int i, j;
    double t;
    NCL*=2;
    
    for (i=0; i<A[0]->m; i++) u[0][i] = 0.0;
    for (i = 0; i < NCL; i+=2) {
        SOR_relax(A[i], f[i], u[i], omega, nr);
        residual(r, f[i], A[i], u[i], 0, 0, NULL);
        dpa_fine2coarse(match[i],   r,      A[i]->n  , f[i+1]); // f[i+1] <- r
        dpa_fine2coarse(match[i+1], f[i+1], A[i+1]->n, f[i+2]); // f[i+2] <- f[i+1]
        for (j = 0; j < A[i+1]->n; j++) u[i+1][j] = 0.0;
        for (j = 0; j < A[i+2]->n; j++) u[i+2][j] = 0.0;
    }
    SOR_relax(A[i], f[i], u[i], omega, nr);
    SOR_relax_rev(A[i], f[i], u[i], omega, nr);
    for (; i > 0; i-=2) {
        dpa_coarse2fine_pure(match[i-1], u[i],   A[i-1]->n, u[i-1]); // u[i-1] <- u[i]
        dpa_coarse2fine     (match[i-2], u[i-1], A[i-2]->n, u[i-2]); // u[i-2] <- u[i-1]
        SOR_relax_rev(A[i-2], f[i-2], u[i-2], omega, nr);
    }
};

void dpa_destroy(int NCL, matrix **A, double **f, double **u, double *r, int **match) {
    int i;
    NCL*=2;
    for (i = 1; i < (NCL + 1); i++) {
        free(f[i]);
        free(u[i]);
    }
    free(f);
    free(u);
    free(r);
    for (i = 1; i < (NCL + 1); i++) destroy_m(A[i]);
    for (i = 1; i < NCL; i++) free(match[i]);
    free(A);
    free(match);
};

int dpa_AMG(matrix *A_o, double *f_o, double *u_o, int NCL, double omega, int nr, double tol, int lmax, double beta) {
    int i, k;
    double **f, **u, *r, delta, t, norm_f;
    matrix **A/*, **I_cf, **I_cf_t*/;
    int **match;

    // AMG_setup gera matrizes A, vetores u, r, e matrizes prolongadoras
    t = dpa_setup(A_o, f_o, u_o, NCL, &A, &f, &u, &r, &match, beta);

    printf("\n\tDPA AMG (%d, %.1f, %d, %.e, %d, %.2f)", NCL, omega, nr, tol, lmax, beta);
    printf("\n\tn/nnz: ");
    for (i = 0; i < ((2 * NCL) + 1); i++)
        if ((i % 2) == 0)
            printf("\n\t\t%d.\t[%d \t%d\t] << SOR", i, A[i]->m, A[i]->nnz);
        else
            printf("\n\t\t%d.\t[%d \t%d\t]", i, A[i]->m, A[i]->nnz);
    printf("\n\tsetup t: %f", t);
    residual(r, f_o, A_o, u_o, 0, 0, NULL);
    norm_f = norm_inf(f_o, A_o->m);
    delta = norm_inf(r, A_o->m) / norm_f;
    k = 0;
    t = get_time();
    while (k < lmax && delta > tol) {
        // multiplica matrizes por u para transferir entre grids
        dpa_Vcycle2(NCL, A, f, u, r, omega, nr, match);
        residual(r, f_o, A_o, u_o, 0, 0, NULL);
        delta = norm_inf(r, A_o->m) / norm_f;
        k++;
        //printf("  %.6f\n", delta);
    }
    t = (get_time() - t) / 100.0;
    dpa_destroy(NCL, A, f, u, r, match);
    return (delta <= tol) ? 1 : 0;
};
