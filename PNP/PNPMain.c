/* PNPMain.c
 * Implementation of enhanced Poisson-Nernst-Plank algorithm proposed by 
 * Dyrka et al (2008)
 * Michelle Shu | July 16, 2013
 */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include "Matrix.h"
# include "SimulationGrid.h"


/* ----- PHYSICAL CONSTANTS ------------------------------------------------ */
const double E    = 1.60217657e-19;  	// Elementary charge (C)
const double E_0  = 8.854187817e-12;  	// Vacuum permittivity (F/m)   
const double K    = 1.3806488e-23;		// Boltzmann constant (J/K)
const float  T    = 293;				// Temperature (K)
const float  E_R  = 78.3;				// Relative permittivity of water

/* ----- CONVERGENCE PARAMETERS -------------------------------------------- */

int main(int argc, char** argv) {
	Simulation* sim = simulation_new();
	print_matrix(sim -> d_na, sim -> M, sim -> N);
	simulation_destroy(sim);
}


/* Update all simulation matrices. Return maximum difference due to update */
double update(Simulation* sim) {
	double max_diff = 0.0;
	// Compute all new ion density distributions
	double* new_ni_na = update_ni(sim -> ni_na, sim -> d_na, sim -> p, 
		sim -> M, sim -> N, 1, sim -> H, &max_diff);
	double* new_ni_k = update_ni(sim -> ni_k, sim -> d_k, sim -> p, 
		sim -> M, sim -> N, 1, sim -> H, &max_diff);
	double* new_ni_cl = update_ni(sim -> ni_cl, sim -> d_cl, sim -> p, 
		sim -> M, sim -> N, -1, sim -> H, &max_diff);
	double* new_ni_a = update_ni(sim -> ni_a, sim -> d_a, sim -> p, 
		sim -> M, sim -> N, -1, sim -> H, &max_diff);
		
	// Compute new electrical potential distribution
	
	return max_diff;	// Return maximum difference
}

/* Compute new ion density matrix based on current values of NI, diffusion
 * coefficients D, and potential P */
double* update_ni(double* ni, double *d, double *p, int M, int N, int Z, 
	int H, double* max_diff) {
	double new_ni[M * N] = {0};
	for (int row = 0; row < N; row++) {
		for (int col = 0; col < M; col++) {
			double numerator = 0.0;
			double denominator = 0.0;
			
			// S and T are intermediate variables to make computations easier
			// Left neighbor
			double S = E * Z / (2 * K * T) * (get_left(p, row, col, M, N) -
				get(p, row, col, M, N));
			double T = (get(d, row, col, M, N) + get_left(d, row, col, M, N)) 
				/ H;
			numerator += get_left(ni) * T * (1 + S);
			denominator += T * (1 - S);
			
			// Right neighbor
			double S = E * Z / (2 * K * T) * (get_right(p, row, col, M, N) -
				get(p, row, col, M, N));
			double T = (get(d, row, col, M, N) + get_right(d, row, col, M, N)) 
				/ H;
			numerator += get_right(ni) * T * (1 + S);
			denominator += T * (1 - S);
			
			// Top neighbor
			double S = E * Z / (2 * K * T) * (get_top(p, row, col, M, N) -
				get(p, row, col, M, N));
			double T = (get(d, row, col, M, N) + get_top(d, row, col, M, N)) 
				/ H;
			numerator += get_top(ni) * T * (1 + S);
			denominator += T * (1 - S);
			
			// Bottom neighbor
			double S = E * Z / (2 * K * T) * (get_bottom(p, row, col, M, N) -
				get(p, row, col, M, N));
			double T = (get(d, row, col, M, N) + get_bottom(d, row, col, M, N)) 
				/ H;
			numerator += get_bottom(ni) * T * (1 + S);
			denominator += T * (1 - S);
			
			// Populate new_ni with new computed value.
			new_ni[get_index(row, col, M)] = numerator / denominator;
			
			// Update max_diff if applicable
			if (abs(get(new_ni, row, col, M, N) - get(ni, row, col, M, N))) >
				*max_diff) {
				*max_diff = abs(get(new_ni, row, col, M, N) - 
					get(ni, row, col, M, N));
			}
		}
	}
}

double* update_potential(Simulation* sim) {
	double* new_p[sim -> M * sim -> N] = {0};
	for (int row = 1; ) 
}