/* Neighbors.h
 * Retrieve desired elements from M x N matrix of doubles internally 
 * represented as 1-D array; indices are ordered in column-major order 
 * (i.e. MATLAB standard)
 *
 * Michelle Shu | July 15, 2013
 */

# include <stdio.h>
# include <stdlib.h>

/* ----- PUBLIC FUNCTION PROTOTYPES ---------------------------------------- */

double get(double *A, int row, int col, int M, int N);
double get_left(double *A, int row, int col, int M, int N);
double get_right(double *A, int row, int col, int M, int N);
double get_top(double *A, int row, int col, int M, int N);
double get_bottom(double *A, int row, int col, int M, int N);
int get_index(int row, int col, int M);
void print_matrix(double *A, int M, int N);
