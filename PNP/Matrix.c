/* Matrix.c
 * Retrieve desired elements from M x N matrix of doubles internally 
 * represented as 1-D array; indices are ordered in column-major order 
 * (i.e. MATLAB standard)
 *
 * Michelle Shu | July 15, 2013
 */

# include "Matrix.h"

/* Test functionality */
/*int main(int argc, char **argv) {
	double A[9] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
	printf("The left neighbor of 5 is %lf\n", getLeft(A, 1, 1, 3, 3));
	printf("The right neighbor of 5 is %lf\n", getRight(A, 1, 1, 3, 3));
	printf("The top neighbor of 5 is %lf\n", getTop(A, 1, 1, 3, 3));
	printf("The bottom neighbor of 5 is %lf\n", getBottom(A, 1, 1, 3, 3));
	printf("The left neighbor of 1 is %lf\n", getLeft(A, 0, 0, 3, 3));
	printf("The right neighbor of 7 is %lf\n", getRight(A, 0, 2, 3, 3));
	printf("The top neighbor of 4 is %lf\n", getTop(A, 0, 1, 3, 3));
	printf("The bottom neighbor of 9 is %lf\n", getBottom(A, 2, 2, 3, 3));
	printf("The right neighbor of 10 is %lf\n", getRight(A, 2, 3, 3, 3));
}*/

/* Get element of matrix A of height M at specified row and column */
double get(double *A, int row, int col, int M, int N) {
	if ((row < 0) || (col < 0) || (row > M - 1) || (col > N - 1)) {
		fprintf(stderr, "Index Out of Bounds: Attempted to retrieve entry "
			"(%d, %d) in %d x %d matrix.\n", row, col, M, N);
		return 0.0;
	}
	return A[(col * M) + row];
}

/* Get left, right, top, bottom neighbors of element (row, col) of matrix A
 * with a height M */
double get_left(double *A, int row, int col, int M, int N) {
	if (col <= 0) {
		/* No wrap, return itself */
		return get(A, row, col, M, N);
	}
	return get(A, row, col - 1, M, N);
}

double get_right(double *A, int row, int col, int M, int N) {
	if (col >= N - 1) {
		/* No wrap, return itself */
		return get(A, row, col, M, N);
	}
	return get(A, row, col + 1, M, N);
}

double get_top(double *A, int row, int col, int M, int N) {
	if (row <= 0) {
		/* No wrap, return itself */
		return get(A, row, col, M, N);
	} else {
		return get(A, row - 1, col, M, N);
	}
}
	
double get_bottom(double *A, int row, int col, int M, int N) {
	if (row >= M - 1) {
		/* No wrap, return itself */
		return get(A, row, col, M, N);
	} else {
		return get(A, row + 1, col, M, N);
	}
}

/* Get array index of element at specified row and column */
int get_index(int row, int col, int M) {
	return (col * M) + row;
}

/* Print matrix contents to console */
void print_matrix(double *A, int M, int N) {
	for (int row = 0; row < M; row++) {
		for (int col = 0; col < N; col++) {
			printf("%e\t", get(A, row, col, M, N));
		}
		printf("\n");
	}
}