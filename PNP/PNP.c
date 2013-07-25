/* PNP.c
 * Implementation of enhanced Poisson-Nernst-Plank algorithm proposed by 
 * Dyrka et al (2008)
 * Michelle Shu | July 16, 2013
 */

# include <float.h>
# include <math.h>
# include <stdbool.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include "Matrix.h"
# include "SimulationGrid.h"
//# include "mex.h"

/* ----- PHYSICAL CONSTANTS ------------------------------------------------ */
const double E    = 1.60217657e-19;  	// Elementary charge (C)
const double E_0  = 8.854187817e-12;  	// Vacuum permittivity (F/m)   
const double K    = 1.3806488e-23;		// Boltzmann constant (J/K)
const float  T    = 293;				// Temperature (K)
const float  E_R  = 78.3;				// Relative permittivity of water

/* ----- CONVERGENCE PARAMETERS -------------------------------------------- */
const float NI_CONV = 1;				// Max allowed ni diff for convergence
const float P_CONV  = 1;				// Max allowed p diff for convergence
const float JUMP_START = 1e-5;
const float JUMP_LOW = 1e-10;			// Lower limit for relaxation step size
const float JUMP_UPP = 1e-2;			// Upper limit for relaxation step size
const int UPDATE_ITER = 100;			// # of iterations to wait per update of
										// jump size

/* ----- PUBLIC FUNCTION PROTOTYPE ----------------------------------------- */

void pnp(Simulation* sim, double* ni_na, double* ni_k, double* ni_cl, 
	double* ni_a, double* ni_p);

/* ----- PRIVATE FUNCTION PROTOTYPES --------------------------------------- */

static void update(Simulation* sim, float jump, double* max_diff_ni, 
		double* max_diff_p);
static void update_ni(Simulation* sim, double* ni, double *d, int Z, 
	double* new_ni, double* max_diff);
static void update_p(Simulation* sim, double* new_p, double* max_diff);

/* ----- PUBLIC FUNCTIONS -------------------------------------------------- */

/*
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	Simulation* sim = simulation_new();
	// Outputs: ni_na, ni_k, ni_cl, ni_a, p 
	plhs[0] = mxCreateDoubleMatrix((mwSize) sim -> M, (mwSize) sim -> N, mxREAL);
	plhs[1] = mxCreateDoubleMatrix((mwSize) sim -> M, (mwSize) sim -> N, mxREAL);
	plhs[2] = mxCreateDoubleMatrix((mwSize) sim -> M, (mwSize) sim -> N, mxREAL);
	plhs[3] = mxCreateDoubleMatrix((mwSize) sim -> M, (mwSize) sim -> N, mxREAL);
	plhs[4] = mxCreateDoubleMatrix((mwSize) sim -> M, (mwSize) sim -> N, mxREAL);
	
	double *ni_na, *ni_k, *ni_cl,*ni_a, *p;
	ni_na = mxGetPr(plhs[0]);
	ni_k = mxGetPr(plhs[1]);
	ni_cl = mxGetPr(plhs[2]);
	ni_a = mxGetPr(plhs[3]);
	p = mxGetPr(plhs[4]);
	
	pnp(sim, ni_na, ni_k, ni_cl, ni_a, p);
	
	memcpy(ni_na, sim -> ni_na, sim -> M * sim -> N * (sizeof(double)));
	memcpy(ni_k, sim -> ni_k, sim -> M * sim -> N * (sizeof(double)));
	memcpy(ni_cl, sim -> ni_cl, sim -> M * sim -> N * (sizeof(double)));
	memcpy(ni_a, sim -> ni_k, sim -> M * sim -> N * (sizeof(double)));
	memcpy(p, sim -> p, sim -> M * sim -> N * (sizeof(double)));

	simulation_destroy(sim);
}

void pnp(Simulation* sim, double* ni_na, double* ni_k, double* ni_cl, 
	double* ni_a, double* p) {
	bool done = false;
	int iter = 0;
	while (! done) {
		done = update(sim);
		if (iter == 1000) {
			done = true;
		}
		iter++;
	}
}

*/
/* C VERSION */
int main(int argc, char** argv) {
	Simulation* sim = simulation_new();
	int iter = 0;
	float jump = JUMP_START;
	double max_diff_ni = NI_CONV;
	double max_diff_p = P_CONV;
	double diff;
	double old_diff = DBL_MAX;	// from last set of iterations

	while (! ((max_diff_ni < NI_CONV) && (max_diff_p < P_CONV))) {
		update(sim, jump, &max_diff_ni, &max_diff_p);
		
		// Adaptive step size adjustment
		// if (iter % UPDATE_ITER == 1) {
		// 	diff = (max_diff_ni + max_diff_p * 100) / 2;
		// 	printf("JUMP: %1.10lf \n", jump);
		// 	printf("MAX_DIFF: %lf \t %lf \n", diff, old_diff);
		// 	printf("MAX_DIFF_P: %lf \n", max_diff_p);
		// 	
		// 	if ((diff > old_diff) && (jump / 2 > JUMP_LOW)) {
		// 		jump /= 2;
		// 	} else if ((diff < old_diff) && (jump * 2 < JUMP_UPP)) {
		// 		jump *= 2;	
		// 	}
		// 	
		// 	old_diff = diff;
		// }
		printf("%lf \t %lf \n", max_diff_ni, max_diff_p);
		if (iter == 10000) {	// Terminate early for testing
			max_diff_ni = 0.0;
			max_diff_p = 0.0;
		}
		iter++;
	}
	for (int col = 0; col < 50; col++) {
		printf("%lf \t %lf \n", get(sim -> p, 20, col, sim -> M, sim -> N), 
		get(sim -> ni_na, 20, col, sim -> M, sim -> N));
	}
	simulation_destroy(sim);
}


/* ----- PRIVATE FUNCTIONS ------------------------------------------------- */

/* Update all simulation matrices. Return true only if max_diff has satisfied
 *convergence conditions */
static void update(Simulation* sim, float jump, double* max_diff_ni, 
	double* max_diff_p) {
		
	double new_ni_na[sim -> M * sim -> N];
	double new_ni_k[sim -> M * sim -> N];
	double new_ni_cl[sim -> M * sim -> N];
	double new_ni_a[sim -> M * sim -> N];
	double new_p[sim -> M * sim -> N];
	
	// Compute all new ion density distributions
	update_ni(sim, sim -> ni_na, sim -> d_na, 1, new_ni_na, max_diff_ni);
	update_ni(sim, sim -> ni_k, sim -> d_k, 1, new_ni_k, max_diff_ni);
	update_ni(sim, sim -> ni_cl, sim -> d_cl, -1, new_ni_cl, max_diff_ni);
	update_ni(sim, sim -> ni_a, sim -> d_a, -1, new_ni_a, max_diff_ni);
		
	// Compute new electrical potential distribution
	update_p(sim, new_p, max_diff_p);
	
	// Relax all values except leftmost and rightmost columns
	for (int row = 0; row < sim -> M; row++) {
		for (int col = 1; col < sim -> N - 1; col++) {
			int i = get_index(row, col, sim -> M);
			sim -> ni_na[i] = ((1 - jump) * sim -> ni_na[i]) + 
				(jump * new_ni_na[i]);
			sim -> ni_k[i] = ((1 - jump) * sim -> ni_k[i]) + 
				(jump * new_ni_k[i]);
			sim -> ni_cl[i] = ((1 - jump) * sim -> ni_cl[i]) + 
				(jump * new_ni_cl[i]);
			sim -> ni_a[i] = ((1 - jump) * sim -> ni_a[i]) + 
				(jump * new_ni_a[i]);
			sim -> p[i] = ((1 - jump) * sim -> p[i]) + (jump * new_p[i]);
		}
	}
}

/* Compute new ion density matrix based on current values of NI, diffusion
 * coefficients D, and potential P */
static void update_ni(Simulation* sim, double* ni, double *d, int Z, 
	double* new_ni, double* max_diff) {
		
	for (int row = 0; row < sim -> M; row++) {
		
		// Static region
		new_ni[get_index(row, 0, sim -> M)] = 0.0;
		new_ni[get_index(row, sim -> N - 1, sim -> M)] = 
			ni[get_index(row, sim -> N - 1, sim -> M)];
		
		// Update region
		for (int col = 1; col < sim -> N - 1; col++) {
			double numerator = 0.0;
			double denominator = 0.0; 
			
			// S and T are intermediate variables to make computations easier
			// Left neighbor
			double S = E * Z / (2 * K * T) * 
				(get_left(sim -> p, row, col, sim -> M, sim -> N) - 
				get(sim -> p, row, col, sim -> M, sim -> N));
			double T = (get(d, row, col, sim -> M, sim -> N) + get_left(d, row, 
				col, sim -> M, sim -> N)) / sim -> H;
			numerator += get_left(ni, row, col, sim -> M, sim -> N) * T * 
				(1 + S);
			denominator += T * (1 - S);
			
			// Right neighbor
			S = E * Z / (2 * K * T) * (get_right(sim -> p, row, col, sim -> M, 
				sim -> N) - get(sim -> p, row, col, sim -> M, sim -> N));
			T = (get(d, row, col, sim -> M, sim -> N) + get_right(d, row, 
				col, sim -> M, sim -> N)) / sim -> H;
			numerator += get_right(ni, row, col, sim -> M, sim -> N) * T * 
				(1 + S);
			denominator += T * (1 - S);

			
			// Top neighbor
			S = E * Z / (2 * K * T) * (get_top(sim -> p, row, col, sim -> M, 
				sim -> N) - get(sim -> p, row, col, sim -> M, sim -> N));
			T = (get(d, row, col, sim -> M, sim -> N) + get_top(d, row, 
				col, sim -> M, sim -> N)) / sim -> H;
			numerator += get_top(ni, row, col, sim -> M, sim -> N) * T * 
				(1 + S);
			denominator += T * (1 - S);

			
			// Bottom neighbor
			S = E * Z / (2 * K * T) * (get_bottom(sim -> p, row, col, sim -> M, 
				sim -> N) - get(sim -> p, row, col, sim -> M, sim -> N));
			T = (get(d, row, col, sim -> M, sim -> N) + get_bottom(d, row, 
				col, sim -> M, sim -> N)) / sim -> H;
			numerator += get_bottom(ni, row, col, sim -> M, sim -> N) * T * 
				(1 + S);
			denominator += T * (1 - S);

			// Set element in new_ni to computed value.
			new_ni[get_index(row, col, sim -> M)] = numerator / denominator;
			
			// Update max_diff if applicable
			if (fabs(get(new_ni, row, col, sim -> M, sim -> N) - 
				get(ni, row, col, sim -> M, sim -> N)) / 
				get(ni, row, col, sim -> M, sim -> N) > *max_diff) {
				*max_diff = fabs(get(new_ni, row, col, sim -> M, sim -> N) - 
					get(ni, row, col, sim -> M, sim -> N)) / 
					get(ni, row, col, sim -> M, sim -> N);
			}
		}
	}
}

static void update_p(Simulation* sim, double* new_p, double* max_diff) {
	for (int row = 0; row < sim -> M; row++) {
		
		// Static region
		new_p[get_index(row, 0, sim -> M)] = 
			sim -> p[get_index(row, 0, sim -> M)];
		new_p[get_index(row, sim -> N - 1, sim -> M)] = 0.0;
		
		// Update region
		for (int col = 1; col < (sim -> N) - 1; col++) {
			// S and T are intermediate variables to make computations easier
			double S = 0.0;
			double T = 0.0;
			
			S += get_left(sim -> p, row, col, sim -> M, sim -> N);	 // Left
			S += get_right(sim -> p, row, col, sim -> M, sim -> N);	 // Right
			S += get_top(sim -> p, row, col, sim -> M, sim -> N);	 // Top
			S += get_bottom(sim -> p, row, col, sim -> M, sim -> N); // Bottom
				
			T += E * get(sim -> ni_na, row, col, sim -> M, sim -> N) / 
				(E_0 * E_R);
			T += E * get(sim -> ni_k, row, col, sim -> M, sim -> N) / 
				(E_0 * E_R);
			T -= E * get(sim -> ni_cl, row, col, sim -> M, sim -> N) / 
				(E_0 * E_R);
			T -= E * get(sim -> ni_a, row, col, sim -> M, sim -> N) / 
				(E_0 * E_R);
			
			// Set element in new_p to computed value.
			new_p[get_index(row, col, sim -> M)] = (S + T) / 4;
			
			// Update max_diff if applicable
			if (fabs(get(new_p, row, col, sim -> M, sim -> N) - 
				get(sim -> p, row, col, sim -> M, sim -> N)) / 
				get(sim -> p, row, col, sim -> M, sim -> N) > *max_diff) {
				*max_diff = fabs(get(new_p, row, col, sim -> M, sim -> N) - 
					get(sim -> p, row, col, sim -> M, sim -> N)) / 
					get(sim -> p, row, col, sim -> M, sim -> N);
			}
		}
	} 
}