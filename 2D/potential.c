/*	POTENTIAL.C
 * 	Evaluate potential distribution within electrical double layer
 *
 *	Michelle Shu | June 28, 2013 
 */

# include <stdio.h>
# include <stdlib.h>
# include <stdbool.h>
# include <math.h>
# include "mex.h"

/* ----- CONSTANTS --------------------------------------------------------- */

const double N_A  = 6.0221413e23;		// Avogadro's number (1/mol)
const double E    = 1.60217657e-19;  	// Elementary charge (C)
const double E_0  = 8.854187817e-12;  	// Vacuum permittivity (F/m)   
const double K    = 1.3806488e-23;		// Boltzmann constant (J/K)
const float  T    = 293;				// Temperature (K)

/* ----- MODEL PARAMETERS -------------------------------------------------- */

const float D    = 2e-8;                // Limit for P to approach 0 
const float L    = 1e-8;                // Length of interface
const float H    = 1e-10;               // Distance step size, determines resolution (m)

const float A	 = 0.1;                 // Potential adjustment ratio
                         	  			// (controls speed of convergence)
const float CONV = 1e-4;             	// Max acceptable error ratio for an individual point for convergence
const float G	 = 1e6;                 // Steepness of initialization curve

/* ----- PUBLIC FUNCTION PROTOTYPE ----------------------------------------- */

void potential(int M, int N, double P_0, double E_R, double EFF, 
			   double ion_types, double *Zi, double *Ci, double MPB, 
			   double *X, double *Y, double *P, double *R);

/* ----- PRIVATE FUNCTION PROTOTYPES --------------------------------------- */

static void initializeGrid(double *X, double *Y, int M, int N, double *P, 
						   double P_0);
static void computePBDerivs(int M, int N, double *P, double *Zi, double *Ci, 
						   	double ion_types, double E_R, double EFF,
							double MPB, double *R);	
static double updatePotentials(int M, int N, double *P, double *R);
static double getNeighborSum(int row, int col, int M, int N, double *P);
static int getIndex(int row, int col, int M);
	
/* ----- FUNCTIONS --------------------------------------------------------- */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	/* Check number of input arguments */
	if (nrhs < 7) {
		mexErrMsgTxt("Not enough input arguments.");
	} else if (nlhs > 4) {
		mexErrMsgTxt("Too many output arguments.");
	}
	
	/* Check type of input arguments */
	for (int i = 0; i < 7; i++) {
		if (!mxIsDouble(prhs[i])) {
			mexErrMsgTxt("Wrong input type.");
		}
	}
	
	int M = (int) (L / H);				// Number of rows
	int N = (int) (D / H) + 1;			// Number of columns
	
	double P_0, E_R, EFF, MPB, ion_types;
	double *Zi, *Ci, *X, *Y, *P, *R;
	
	/* Parse input arguments */
	P_0 = mxGetScalar(prhs[0]);
	E_R = mxGetScalar(prhs[1]);
	EFF = mxGetScalar(prhs[2]);
	ion_types = mxGetScalar(prhs[3]);
	Zi = mxGetPr(prhs[4]);
	Ci = mxGetPr(prhs[5]);
	MPB = mxGetScalar(prhs[6]);
	
	/* Set output arguments */
	plhs[0] = mxCreateDoubleMatrix((mwSize) 1, (mwSize) N, mxREAL); // X
	plhs[1] = mxCreateDoubleMatrix((mwSize) 1, (mwSize) M, mxREAL); // Y
	plhs[2] = mxCreateDoubleMatrix((mwSize) M, (mwSize) N, mxREAL); // P
	plhs[3] = mxCreateDoubleMatrix((mwSize) M, (mwSize) N, mxREAL); // R
	
	X = mxGetPr(plhs[0]);
	Y = mxGetPr(plhs[1]);
	P = mxGetPr(plhs[2]);
	R = mxGetPr(plhs[3]);
	
	potential(M, N, P_0, E_R, EFF, ion_types, Zi, Ci, MPB, X, Y, P, R);
}

void potential(int M, int N, double P_0, double E_R, double EFF, 
			   double ion_types, double *Zi, double *Ci, double MPB, 
			   double *X, double *Y, double *P, double *R) {
	
	initializeGrid(X, Y, M, N, P, P_0);
	
	double max_error = CONV + 1;
	while (max_error > CONV) {
		computePBDerivs(M, N, P, Zi, Ci, ion_types, E_R, EFF, MPB, R);
		max_error = updatePotentials(M, N, P, R);
	}
}


static void initializeGrid(double *X, double *Y, int M, int N, double *P, 
						   double P_0) {	
	
	double axis_val = 0.0;
	for (int i = 0; i < N; i++) {
		X[i] = axis_val;
		axis_val += H;
	}
	
	axis_val = 0.0;
	for (int i = 0; i < M; i++) {
		Y[i] = axis_val;
		axis_val += H;
	}
	
	for (int col = 0; col < N; col++) {
		double potential = (- P_0 / D * X[col]) + P_0;
		// double potential = (P_0 / (1 - exp(-G * D))) * 
		//				   (exp(-G * X[col]) - exp(-G * D));
		for (int row = 0; row < M; row++) {
			P[getIndex(row, col, M)] = potential;
		}
	}
}

/** Use (Modified) Poisson Boltzmann equation to compute second derivatives for
 *  each element of P */
static void computePBDerivs(int M, int N, double *P, double *Zi, double *Ci, 
							double ion_types, double E_R, double EFF, 
							double MPB, double *R) {
								
	for (int row = 0; row < M; row++) {
		for (int col = 0; col < N; col++){
			int i = getIndex(row, col, M);
			R[i] = 0;
			for (int ion_num = 0; ion_num < ion_types; ion_num++) {
				double Z = Zi[ion_num];
				double C_0 = Ci[ion_num];
				double V = 2 * (pow(EFF, 3)) * C_0 * N_A;
				
				double s1 = sinh(Z * E * P[i] / (K * T));
				R[i] += (2 * Z * E * N_A * C_0 / (E_0 * E_R)) * s1;
				if (MPB == 1) {
					double s2 = sinh(Z * E * P[i] / (2 * K * T));
					R[i] /= 1 + (2 * V * s2 * s2);
				}
			}
		}
	}
}

/** Update potential values based on PB second derivative and potentials of
 * neighbors. Return the maximum deviation (error) encountered. */
static double updatePotentials(int M, int N, double *P, double *R) {
	
	double max_error = 0;
	for (int row = 0; row < M; row++) {
		for (int col = 1; col < N - 1; col++) {	// Fix first and last points
			int i = getIndex(row, col, M);
			double new_p = (getNeighborSum(row, col, M, N, P) - 
						   (R[i] * H * H)) / 4;
			// Evaluate error
			double error = fabs(new_p - P[i]) / P[i];
			if (error > max_error) {
				max_error = error;
			}
			
			// Update P
			P[i] += A * (new_p - P[i]);
		}
	}
	
	return max_error;
}

/** Get the sum of the potentials of all immediate neighboring locations */
static double getNeighborSum(int row, int col, int M, int N, double *P) {
	
	double neighbor_sum = 0.0;
	
	// Left neighbor
	if (col > 0) {
		neighbor_sum += P[getIndex(row, col - 1, M)];
	} else {
		neighbor_sum += P[getIndex(row, N - 1, M)];
	}
	
	// Right neighbor
	if (col < N - 1) {
		neighbor_sum += P[getIndex(row, col + 1, M)];
	} else {
		neighbor_sum += P[getIndex(row, 0, M)];
	}
	
	// Top neighbor
	if (row > 0) {
		neighbor_sum += P[getIndex(row - 1, col, M)];
	} else {
		neighbor_sum += P[getIndex(M - 1, col, M)];
	}
	
	// Bottom neighbor
	if (row < M - 1) {
		neighbor_sum += P[getIndex(row + 1, col, M)];
	} else {
		neighbor_sum += P[getIndex(0, col, M)];
	}
	
	return neighbor_sum;
}

static int getIndex(int row, int col, int M) {
	return (col * M) + row;
}