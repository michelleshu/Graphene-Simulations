/*	potential.c
 * 	Evaluate potential distribution within electrical double layer
 *
 *	Michelle Shu | June 23, 2013 
 */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

/* ----- CONSTANTS ----------------------------------------------------------------- */

const double N_A  = 6.0221413e23;		// Avogadro's number (1/mol)
const double E    = 1.60217657e-19;  	// Elementary charge (C)
const double E_0  = 8.854187817e-12;  	// Vacuum permittivity (F/m)   
const double K    = 1.3806488e-23;		// Boltzmann constant (J/K)
const float T     = 293;				// Temperature (K)

/* ----- MODEL PARAMETERS ---------------------------------------------------------- */

const float D    = 3e-8;                // Limit for P to approach 0 (Max X value)
const float L    = 1e-8;                // Length of interface
const float H    = 1e-10;               // Distance step size, determines resolution (m)

const float A	  = 0.1;                 // Initial potential function adjustment step ratio
                         	  			 // (controls speed of convergence)
const float CONV = 1e-4;             	 // Convergence criterion: Max acceptable ratio of
                         	   			 // (P_calc - P)/ P for an individual point
const float G	  = 1e6;                 // Steepness of initialization curve

/* ----- PRIVATE FUNCTION PROTOTYPES ----------------------------------------------- */

static void initializeGrid(double *X, double *Y, int M, int N, double *P, float P_0);
static double getValueAt(int row, int col, int N, double *P);
static int getIndexFor(int row, int col, int N);
static int getRowOf(int index, int N);
static int getColOf(int index, int N);
	
/* ----- FUNCTIONS ----------------------------------------------------------------- */

int main(int argc, char *argv[]) {
	int M = (int) (L / H);			// Number of rows
	int N = (int) (D / H) + 1;		// Number of columns
	
	float P_0 = 0.025;
	int Z = 1;
	int C = 150;
	float E_R = 78.3;
	float EFF = 1e-9;
	int MPB = 1;

	double X[N];					// Array of x values
	double Y[M];					// Array of y values
	double P[M * N];				// Array of potential values
	double R[M * N];				// Array of potential second derivatives based on PB eq.
	double P_calc[M * N] = {0};		// Array of adjusted P values
	
	initializeGrid(X, Y, M, N, P, P_calc, P_0);
	
	double max_error = CONV + 1;
	while (max_error > CONV) {
		computePBDerivs
	}
}


static void initializeGrid(double *X, double *Y, int M, int N, double *P, double *P_calc, float P_0) {	
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
		double potential = (P_0 / (1 - exp(-G * D))) * (exp(-G * X[col]) - exp(-G * D));
		for (int row = 0; row < M; row++) {
			P[getIndexFor(row, col, N)] = potential;
		}
	}
	
	for (int row = 0; row < M; row++) {
		P_calc[getIndexFor(row, 0, N)] = P_0;
		P_calc[getIndexFor(row, M - 1, N)] = 0;
	}
}

/** Use (Modified) Poisson Boltzmann equation to compute second derivatives for each element of P */
static void computePBDerivs(int M, int N, double *P, int Z, int C, float E_R, float EFF, int MPB) {
	
}

static double getValueAt(int row, int col, int N, double *A) {
	return A[(row * N) + col];
}

static int getIndexFor(int row, int col, int N) {
	return (row * N) + col;
}

static int getRowOf(int index, int N) {
	return (int) (index / N);
}

static int getColOf(int index, int N) {
	return (index % N);
}