# ifndef bcm_h
# define bcm_h

#include "matrix.h"

int bcm_trymatch(int rowindex, int colindex, matrix *W, int *jrowindex,
        int ljrowindex, int *jcolindex, int ljcolindex, int *rmatch);
int * bcm_CSRMatrixHMatch(matrix *A);

# endif
