/* Some definitions and notes:
Matrix Data Structure: 	A continuous single dimension array of concatenated rows in order represented as row 1 + row 2 + row 3 + ....
Augmented Matrix: 		A matrix in the form of linear equations, with the variables(coefficients) in the left columns and the solutions(constants) in the rightmost last column.
*/


#include <stdio.h>
#include <stdlib.h>

#include "matrixalgebra.h"


/* Matrix Util Functions */

// Matrix and Matrix->array must be freed (matrix_free does both)
Matrix* matrix_create(M_TYPE* matrixArray, size_t rows, size_t columns) {
	Matrix *newMatrix = (Matrix *) malloc(sizeof(Matrix));
	newMatrix->row = rows;
	newMatrix->col = columns;

	size_t arraySize = rows * columns;
	newMatrix->array = (M_TYPE *) malloc(sizeof(M_TYPE) * arraySize);
	for (size_t i = 0; i < arraySize; ++i) {
		(newMatrix->array)[i] = matrixArray[i];
	}

	return newMatrix;
}

// Create diagonal square matrix of dimSize x dimSize with passed valued. Must be free'ed with matrix_free.
Matrix* matrix_create_diagonal(M_TYPE diagonalValue, size_t dimSize) {
	Matrix *newMatrix = (Matrix *) malloc(sizeof(Matrix));
	newMatrix->row = dimSize;
	newMatrix->col = dimSize;

	newMatrix->array = (M_TYPE *) calloc(dimSize * dimSize, sizeof(M_TYPE));
	for (size_t i = 0; i < dimSize; ++i) {
		(newMatrix->array)[(i * dimSize) + i] = diagonalValue;
	}

	return newMatrix;
}

void matrix_free(Matrix* matrix) {
	free(matrix->array);
	free(matrix);
}

void matrix_print(Matrix* matrix) {
	printf("Matrix Size: R:%zu | C:%zu\n", matrix->row, matrix->col);
	for(size_t i = 0; i < matrix->row; ++i) {
		for(size_t j = 0; j < matrix->col; ++j) {
			#if __STDC_VERSION__ >= 201112L
				printf("\t");
				printNumber(matrix->array[(i * matrix->col) + j]);
			#else
				printf("\t%Lf", matrix->array[(i * matrix->col) + j]);
			#endif
		}
		printf("\n");
	}
}


/* Elementary Row Operations */

// Scalar Multiplication: s * R1 -> R1
void eroScale(Matrix* matrix, size_t rowIndex, S_TYPE scalar) {
	size_t rowOffset = rowIndex * matrix->col;
	for (size_t columnOffset = 0; columnOffset < matrix->col; ++columnOffset) {

		if (!ENTRY_COMPARE((matrix->array)[rowOffset + columnOffset], (M_TYPE) 0)) {
			(matrix->array)[rowOffset + columnOffset] = (M_TYPE) (((double) ((matrix->array)[rowOffset + columnOffset])) * scalar);
		}
	}
}

// Swap Interchange: R1 <-> R2
void eroSwap(Matrix* matrix, size_t rowIndex1, size_t rowIndex2) {
	size_t rowOffset1 = rowIndex1 * matrix->col;
	size_t rowOffset2 = rowIndex2 * matrix->col;
	for (size_t columnOffset = 0; columnOffset < matrix->col; ++columnOffset) {
		M_TYPE temp = (matrix->array)[rowOffset1 + columnOffset];
		(matrix->array)[rowOffset1 + columnOffset] = (matrix->array)[rowOffset2 + columnOffset];
		(matrix->array)[rowOffset2 + columnOffset] = temp;
	}
}

// Scale Sum: sR1 + R2 -> R2
void eroSum(Matrix* matrix, size_t rowIndex1, size_t rowIndex2, S_TYPE scalar) {
	size_t rowOffset1 = rowIndex1 * matrix->col;
	size_t rowOffset2 = rowIndex2 * matrix->col;
	for (size_t columnOffset = 0; columnOffset < matrix->col; ++columnOffset) {
		(matrix->array)[rowOffset2 + columnOffset] += (M_TYPE) (((S_TYPE) ((matrix->array)[rowOffset1 + columnOffset])) * scalar);
	}
}


/* Echelon Form */

// REF: reducedRowEchelonForm = 0, RREF: reducedRowEchelonForm = 1
void rowEchelonForm(Matrix* matrix, int reducedRowEchelonForm) {
	size_t currentWorkingRow = 0; // Track which row is next to have leading one

	for (size_t column = 0; column < matrix->col - 1; ++column) { // For every column except last (solution) column

		size_t lastNonZeroRowIndex = 0;
		int isZeroColumn = 1; // Assume column is filled with zeros

		size_t leadingRowOffset = currentWorkingRow * matrix->col;
		size_t leadingRowCoefficientOffset = leadingRowOffset + column;

		if (ENTRY_COMPARE((matrix->array)[leadingRowCoefficientOffset], (M_TYPE) 0)) { // If leading coefficient zero, swap for non-zero or skip column

			if (currentWorkingRow < matrix->row - 1) { // Avoid looping when on last row
				for (size_t row = currentWorkingRow; row < matrix->row; ++row) { // Find first non-zero row and swap it
					size_t rowOffset = row * matrix->col;
					size_t totalOffset = rowOffset + column;
					if (!ENTRY_COMPARE((matrix->array)[totalOffset], (M_TYPE) 0)) { // If non-zero found
						isZeroColumn = 0; // No longer zero column
						lastNonZeroRowIndex = row;
						eroSwap(matrix, currentWorkingRow, row); // Swap with non-zero row
						break;
					}
				}
			}

			if (isZeroColumn) { // If all rows have zeros, ignore this column
				continue;
			}
		}

		if ((matrix->array)[leadingRowCoefficientOffset] != (M_TYPE) 1) { // If leading coefficient not one, scale it to one
			eroScale(matrix, currentWorkingRow, ((S_TYPE) 1)/((S_TYPE) (matrix->array)[leadingRowCoefficientOffset]));
		}

		// currentWorkingRow now has leading one

		// Bottom rows to zeros
		size_t startingBottomRow = (isZeroColumn) ? currentWorkingRow : lastNonZeroRowIndex;
		for (size_t row = startingBottomRow + 1; row < matrix->row; ++row) { // Loop to convert bottom to zero
			if (!ENTRY_COMPARE((matrix->array)[(row * matrix->col) + column], (M_TYPE) 0)) {
				eroSum(matrix, currentWorkingRow, row, ((S_TYPE) -1) * ((S_TYPE) (matrix->array)[(row * matrix->col) + column]));
			}
		}

		// Top rows to zeros
		if (reducedRowEchelonForm && currentWorkingRow) {
			for (size_t offsetRow = currentWorkingRow; offsetRow > 0; --offsetRow) {
				size_t row = offsetRow - 1;
				if (!ENTRY_COMPARE((matrix->array)[(row * matrix->col) + column], (M_TYPE) 0)) {
					eroSum(matrix, currentWorkingRow, row, ((S_TYPE) -1) * ((S_TYPE) (matrix->array)[(row * matrix->col) + column]));
				}
			}
		}

		++currentWorkingRow;
	}
}

// Convert augmented matrix to Row Echelon Form
void matrix_ref(Matrix* matrix) {
	rowEchelonForm(matrix, 0);
}

// Convert augmented matrix to Reduced Row Echelon Form
void matrix_rref(Matrix* matrix) {
	rowEchelonForm(matrix, 1);
}


/* Existence and Uniqueness of Solutions */

// Check if rref augmented matrix has no solutions (UOS_ZERO), one solution (UOS_ONE), or infinite solutions (UOS_INF)
// Assumes all zero rows are at the bottom
int matrix_uniquenessOfSolutions(Matrix* matrix) {
	int isConsistent = 0; // Assume bottom zero row exists

	// Check if inconsistent
	for (size_t offsetRow = matrix->row; offsetRow > 0; --offsetRow) {

		size_t row = offsetRow - 1;

		if (!ENTRY_COMPARE((matrix->array)[((row + 1) * matrix->col) - 1], (M_TYPE) 0)) {

			for (size_t offsetColumn = matrix->col - 1; offsetColumn > 0; --offsetColumn) {
				size_t column = offsetColumn - 1; // (0 start index + augmented solution column at end) compensation
				if (!ENTRY_COMPARE((matrix->array)[(row * matrix->col) + column], (M_TYPE) 0)) {
					isConsistent = 1;
					break;
				}
			}

			if (!isConsistent) { // If row is all zeros, return inconsistent solution
				return UOS_ZERO;
			}

			break;
		}
	}


	if (!isConsistent) { // Matrix is all zeros
		return UOS_INF;
	}

	if (matrix->row < (matrix->col - 1)) { // More variables than equations (but could be that some equations are all zeros)
		return UOS_INF;
	}

	// Check if consistent and unique
	for (size_t column = 0; column < (matrix->col - 1); ++column) {

		int leadingOne = 0;

		for (size_t row = 0; row < matrix->row; ++row) {
			if (!ENTRY_COMPARE((matrix->array)[(row * matrix->col) + column], (M_TYPE) 0)) {
				if (leadingOne) {
					return UOS_INF; // Ruh roh Raggy, there be more than one leading one.
				}
				leadingOne = 1;
			}
		}
	}

	return UOS_ONE;
}


/* Matrix Operations */

// Returned matrices must be free'ed, assuming they aren't NULL

// Returns 1 if equal dimensions, 0 otherwise
int matrix_isEqualDimensions(Matrix* matrix1, Matrix* matrix2) {
	if ((matrix1->row != matrix2->row) || (matrix1->col != matrix2->col)) {
		return 0;
	}

	return 1;
}

// Returns 1 if matrices are equivalent, 0 otherwise
int matrix_isEqual(Matrix* matrix1, Matrix* matrix2) {
	if (matrix_isEqualDimensions(matrix1, matrix2)) {
		return 0;
	}

	for (size_t arrayIndex = 0; arrayIndex < ((matrix1->row * matrix1->col) - 1); ++arrayIndex) {
		if ((matrix1->array)[arrayIndex] != (matrix2->array)[arrayIndex]) {
			return 0;
		}
	}

	return 1;
}

// Return 1 if zero matrix, 0 otherwise
int matrix_isZero(Matrix* matrix) {
	for (size_t arrayIndex = 0; arrayIndex < ((matrix->row * matrix->col) - 1); ++arrayIndex) {
		if (!ENTRY_COMPARE((matrix->array)[arrayIndex], (M_TYPE) 0)) {
			return 0;
		}
	}

	return 1;
}

// Matrix1 + Matrix2
Matrix* matrix_add(Matrix* matrix1, Matrix* matrix2) {
	if (!matrix_isEqualDimensions(matrix1, matrix2)) {
		return NULL;
	}

	size_t matrixSize = matrix1->row * matrix1->col;

	M_TYPE *sumMatrixVector = (M_TYPE *) malloc(sizeof(M_TYPE) * matrixSize);

	for (size_t arrayIndex = 0; arrayIndex < (matrixSize - 1); ++arrayIndex) {
		sumMatrixVector[arrayIndex] = (matrix1->array)[arrayIndex] + (matrix2->array)[arrayIndex];
	}

	Matrix *sumMatrix = matrix_create(sumMatrixVector, matrix1->row, matrix1->col);
	free(sumMatrixVector);

	return sumMatrix;
}

// Scalar * Matrix
Matrix* matrix_scale(Matrix* matrix, S_TYPE scalar) {
	Matrix *scaledMatrix = matrix_create(matrix->array, matrix->row, matrix->col);

	for (size_t row = 0; row < (scaledMatrix->row - 1); ++row) {
		eroScale(scaledMatrix, row, scalar);
	}

	return scaledMatrix;
}

// Matrix1 * Matrix2
Matrix* matrix_compose(Matrix* matrix1, Matrix* matrix2) {
	if (matrix1->col != matrix2->row) { // Inner dimension check
		return NULL;
	}

	M_TYPE *productMatrixVector = (M_TYPE *) calloc(matrix1->row * matrix2->col, sizeof(M_TYPE));

	for (size_t column2Index = 0; column2Index < matrix2->col; ++column2Index) { // For every matrix2 column
		for (size_t row2Index = 0; row2Index < matrix2->row; ++row2Index) { // Go through each scalar

			M_TYPE scalar = (matrix2->array)[(row2Index * matrix2->col) + column2Index];

			for (size_t row1Index = 0; row1Index < matrix1->row; ++row1Index) {

				M_TYPE m1Entry = (matrix1->array)[(row1Index * matrix1->col) + row2Index];

				productMatrixVector[(row1Index * matrix2->col) + column2Index] += scalar * m1Entry;
			}

		}
	}

	Matrix *productMatrix = matrix_create(productMatrixVector, matrix1->row, matrix2->col);
	free(productMatrixVector);

	return productMatrix;
}

// [Matrix1 Matrix2] augmented matrix
Matrix* matrix_augment(Matrix* matrix1, Matrix* matrix2) {
	if (matrix1->row != matrix2->row) {
		return NULL;
	}

	Matrix *augmentedMatrix = (Matrix *) malloc(sizeof(Matrix));
	augmentedMatrix->row = matrix1->row;
	augmentedMatrix->col = matrix1->col + matrix2->col;

	augmentedMatrix->array = (M_TYPE *) malloc(sizeof(M_TYPE) * augmentedMatrix->row * augmentedMatrix->col);

	for (size_t row = 0; row < matrix1->row; ++row) {

		for (size_t col = 0; col < matrix1->col; ++col) {
			(augmentedMatrix->array)[(row * (matrix1->col + matrix2->col)) + col] = (matrix1->array)[(row * matrix1->col) + col];
		}

		for (size_t col = 0; col < matrix2->col; ++col) {
			(augmentedMatrix->array)[(row * (matrix1->col + matrix2->col)) + (matrix1->col + col)] = (matrix2->array)[(row * matrix2->col) + col];
		}

	}

	return augmentedMatrix;
}

// AX = B, return X matrix solution, or NULL if inf/none
Matrix* matrix_solve(Matrix* matrixA, Matrix* matrixB) {
	if (matrixA->row != matrixA->col) {
		return NULL;
	}

	Matrix* augmentedMatrix = matrix_augment(matrixA, matrixB);

	if (!augmentedMatrix) {
		return NULL;
	}

	matrix_rref(augmentedMatrix);

	// Check if rref form has identity matrix in place of matrixA
	for (size_t rowA = 0; rowA < matrixA->row; ++rowA) {
		if (!ENTRY_COMPARE((augmentedMatrix->array)[rowA * (augmentedMatrix->col + 1)], (M_TYPE) 1)) {
			return NULL;
		}
	}

	Matrix* xMatrix = (Matrix *) malloc(sizeof(Matrix));
	xMatrix->row = matrixB->row;
	xMatrix->col = matrixB->col;
	xMatrix->array = (M_TYPE *) malloc(sizeof(M_TYPE) * xMatrix->row * xMatrix->col);

	// Copy X matrix
	for (size_t rowB = 0; rowB < matrixB->row; ++rowB) {
		for (size_t colB = 0; colB < matrixB->col; ++colB) {
			(xMatrix->array)[(rowB * xMatrix->col) + colB] = (augmentedMatrix->array)[(rowB * augmentedMatrix->col) + (colB + matrixA->col)];
		}
	}

	matrix_free(augmentedMatrix);

	return xMatrix;
}

Matrix* matrix_transpose(Matrix* matrix) {
	Matrix* transposedMatrix = (Matrix *) malloc(sizeof(Matrix));
	transposedMatrix->row = matrix->col;
	transposedMatrix->col = matrix->row;
	transposedMatrix->array = (M_TYPE *) malloc(sizeof(M_TYPE) * matrix->row * matrix->col);

	for (size_t row = 0; row < matrix->row; ++row) {
		for (size_t col = 0; col < matrix->col; ++col) {
			(transposedMatrix->array)[(col * transposedMatrix->col) + row] = (matrix->array)[(row * matrix->col) + col];
		}
	}

	return transposedMatrix;
}

// Returns 0 is non-square matrix
M_TYPE matrix_trace(Matrix* matrix) {
	if (matrix->row != matrix->col) {
		return 0;
	}

	M_TYPE trace = 0;
	for (size_t i = 0; i < matrix->row; ++i) {
		trace += (matrix->array)[(i * matrix->col) + i];
	}

	return trace;
}

// Returns 0 is non-square matrix
M_TYPE matrix_determinant(Matrix* matrix) {
	if (matrix->row != matrix->col) {
		return 0;
	}

	if (matrix->row == 1) {
		return (matrix->array)[0];
	}

	if (matrix->row == 2) {
		M_TYPE a = (matrix->array)[0];
		M_TYPE b = (matrix->array)[1];
		M_TYPE c = (matrix->array)[2];
		M_TYPE d = (matrix->array)[3];

		return ((a * d) - (b * c));
	}

	M_TYPE det = (M_TYPE) 0;
	for (size_t j = 0; j < matrix->col; ++j) {
		M_TYPE scalar = (matrix->array)[j];

		Matrix *minor = (Matrix *) malloc(sizeof(Matrix));
		minor->row = matrix->row - 1;
		minor->col = matrix->col - 1;
		minor->array = (M_TYPE *) malloc(sizeof(M_TYPE) * minor->row * minor->col);

		int columnOffset = 0;

		// Copy minor of matrix to new matrix
		for (size_t column = 0; column < matrix->col; ++column) {
			// Skip if it's the current column j
			if (column == j) {
				columnOffset++;
				continue;
			}

			for (size_t row = 1; row < matrix->row; ++row) {
				(minor->array)[((row - 1) * minor->col) + (column - columnOffset)] = (matrix->array)[(row * matrix->col) + column];
			}
		}

		M_TYPE matrixMinor = matrix_determinant(minor);
		matrix_free(minor);

		if (j & 1) { // Odd exponent (i+j, which is really 1+j, but j starts at 0, so it's really 1 + j + 1, which can just be j)
			det -= (matrixMinor * scalar);
		} else { // Even
			det += (matrixMinor * scalar);
		}


	}

	return det;
}