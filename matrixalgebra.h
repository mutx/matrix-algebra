#define M_TYPE float // Matrix entry type
#define S_TYPE double // Scalar type

/* Uniqueness of Solutions */
#define UOS_ZERO	0
#define UOS_ONE		1
#define UOS_INF		2

/* Compare two matrice values */
#define ENTRY_COMPARE(x, y) ((x) == (y) ? 1 : 0) // Requires epsilon approximation fix

/* http://www.robertgamble.net/2012/01/c11-generic-selections.html */
#if __STDC_VERSION__ >= 201112L
	#define printf_dec_format(x) _Generic((x), \
	    char: "%c", \
	    signed char: "%hhd", \
	    unsigned char: "%hhu", \
	    signed short: "%hd", \
	    unsigned short: "%hu", \
	    signed int: "%d", \
	    unsigned int: "%u", \
	    long int: "%ld", \
	    unsigned long int: "%lu", \
	    long long int: "%lld", \
	    unsigned long long int: "%llu", \
	    float: "%f", \
	    double: "%f", \
	    long double: "%Lf", \
	    char *: "%s", \
	    void *: "%p")

	#define printNumber(x) printf(printf_dec_format(x), x);
#endif

typedef struct matrix {
	size_t row, col; // Row, Column
	M_TYPE *array; // Pointer to preallocated array filled with values
} Matrix;


Matrix* matrix_create				(M_TYPE* matrixArray, size_t rows, size_t columns);
Matrix* matrix_create_diagonal		(M_TYPE diagonalValue, size_t dimSize);
void 	matrix_free					(Matrix* matrix);
void 	matrix_print				(Matrix* matrix);

void 	eroScale					(Matrix* matrix, size_t rowIndex, S_TYPE scalar);
void 	eroSwap						(Matrix* matrix, size_t rowIndex1, size_t rowIndex2);
void 	eroSum						(Matrix* matrix, size_t rowIndex1, size_t rowIndex2, S_TYPE scalar);
void 	rowEchelonForm				(Matrix* matrix, int reducedRowEchelonForm);
void 	matrix_ref					(Matrix* matrix);
void 	matrix_rref					(Matrix* matrix);

int 	matrix_uniquenessOfSolutions(Matrix* matrix);
int 	matrix_isEqualDimensions	(Matrix* matrix1, Matrix* matrix2);
int 	matrix_isZero				(Matrix* matrix);

Matrix* matrix_add					(Matrix* matrix1, Matrix* matrix2);
Matrix* matrix_scale				(Matrix* matrix, S_TYPE scalar);
Matrix* matrix_compose				(Matrix* matrix1, Matrix* matrix2);
Matrix* matrix_augment				(Matrix* matrix1, Matrix* matrix2);
Matrix* matrix_solve				(Matrix* matrixA, Matrix* matrixB);
Matrix* matrix_transpose			(Matrix* matrix);
M_TYPE 	matrix_trace				(Matrix* matrix);
M_TYPE 	matrix_determinant			(Matrix* matrix);