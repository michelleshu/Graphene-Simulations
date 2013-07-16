/* PNP.c
 * Implementation of enhanced Poisson-Nernst-Plank algorithm proposed by 
 * Dyrka et al (2008)
 * Michelle Shu | July 16, 2013
 */

# include <math.h>
# include <stdbool.h>
# include <stdio.h>
# include <stdlib.h>
# include "Matrix.h"
# include "SimulationGrid.h"

/* ----- PHYSICAL CONSTANTS ------------------------------------------------ */
const double E    = 1.60217657e-19;  	// Elementary charge (C)
const double E_0  = 8.854187817e-12;  	// Vacuum permittivity (F/m)   
const double K    = 1.3806488e-23;		// Boltzmann constant (J/K)
const float  T    = 293;				// Temperature (K)
const float  E_R  = 78.3;				// Relative permittivity of water

/* ----- CONVERGENCE PARAMETERS -------------------------------------------- */
const float NI_CONV = 0.01;			// Max allowed ni diff for convergence
const float P_CONV  = 0.01;			// Max allowed p diff for convergence
const float JUMP	= 0.01;				// Jump size for relaxation step

/* ----- PRIVATE FUNCTION PROTOTYPES --------------------------------------- */

static bool update(Simulation* sim);
static double* update_ni(double* ni, double *d, double *p, int M, int N, int Z, 
	int H, double* max_diff);
static double* update_p(Simulation* sim, double* max_diff);

/* ----- PUBLIC FUNCTIONS -------------------------------------------------- */

int main(int argc, char** argv) {
	Simulation* sim = simulation_new();
	print_matrix(sim -> d_na, sim -> M, sim -> N);
	
	bool done = false;
	while (! done) {
		done = update(sim);
	}
	print_matrix(sim -> d_na, sim -> M, sim -> N);
	
	simulation_destroy(sim);
}

/* ----- PRIVATE FUNCTIONS ------------------------------------------------- */

/* Update all simulation matrices. Return true only if max_diff has satisfied
 *convergence conditions */
static bool update(Simulation* sim) {
	double* max_diff_ni;
	double* max_diff_p;
	
	// Compute all new ion density distributions
	double* new_ni_na = update_ni(sim -> ni_na, sim -> d_na, sim -> p, 
		sim -> M, sim -> N, 1, sim -> H, max_diff_ni);
	double* new_ni_k = update_ni(sim -> ni_k, sim -> d_k, sim -> p, 
		sim -> M, sim -> N, 1, sim -> H, max_diff_ni);
	double* new_ni_cl = update_ni(sim -> ni_cl, sim -> d_cl, sim -> p, 
		sim -> M, sim -> N, -1, sim -> H, max_diff_ni);
	double* new_ni_a = update_ni(sim -> ni_a, sim -> d_a, sim -> p, 
		sim -> M, sim -> N, -1, sim -> H, max_diff_ni);
		
	// Compute new electrical potential distribution
	double* new_p = update_p(sim, max_diff_p);
	
	// Relax (update) all values
	for (int i = 0; i < (sim -> M * sim -> N); i++) {
		sim -> ni_na[i] = ((1 - JUMP) * sim -> ni_na[i]) + (JUMP * new_ni_na[i]);
		sim -> ni_k[i] = ((1 - JUMP) * sim -> ni_k[i]) + (JUMP * new_ni_k[i]);
		sim -> ni_cl[i] = ((1 - JUMP) * sim -> ni_cl[i]) + (JUMP * new_ni_cl[i]);
		sim -> ni_a[i] = ((1 - JUMP) * sim -> ni_a[i]) + (JUMP * new_ni_a[i]);
		sim -> p[i] = ((1 - JUMP) * sim -> p[i]) + (JUMP * new_p[i]);
	}
	
	return ((*max_diff_ni <= NI_CONV) && (*max_diff_p <= P_CONV));
}

/* Compute new ion density matrix based on current values of NI, diffusion
 * coefficients D, and potential P */
static double* update_ni(double* ni, double *d, double *p, int M, int N, int Z, 
	int H, double* max_diff) {
	double new_ni[M * N];
	
	for (int row = 0; row < M; row++) {
		for (int col = 0; col < N; col++) {
			double numerator = 0.0;
			double denominator = 0.0;
			
			// S and T are intermediate variables to make computations easier
			// Left neighbor
			double S = E * Z / (2 * K * T) * (get_left(p, row, col, M, N) -
				get(p, row, col, M, N));
			double T = (get(d, row, col, M, N) + get_left(d, row, col, M, N)) 
				/ H;
			numerator += get_left(ni, row, col, M, N) * T * (1 + S);
			denominator += T * (1 - S);
			
			// Right neighbor
			S = E * Z / (2 * K * T) * (get_right(p, row, col, M, N) -
				get(p, row, col, M, N));
			T = (get(d, row, col, M, N) + get_right(d, row, col, M, N)) 
				/ H;
			numerator += get_right(ni, row, col, M, N) * T * (1 + S);
			denominator += T * (1 - S);
			
			// Top neighbor
			S = E * Z / (2 * K * T) * (get_top(p, row, col, M, N) -
				get(p, row, col, M, N));
			T = (get(d, row, col, M, N) + get_top(d, row, col, M, N)) 
				/ H;
			numerator += get_top(ni, row, col, M, N) * T * (1 + S);
			denominator += T * (1 - S);
			
			// Bottom neighbor
			S = E * Z / (2 * K * T) * (get_bottom(p, row, col, M, N) -
				get(p, row, col, M, N));
			T = (get(d, row, col, M, N) + get_bottom(d, row, col, M, N)) 
				/ H;
			numerator += get_bottom(ni, row, col, M, N) * T * (1 + S);
			denominator += T * (1 - S);
			
			// Set element in new_ni to computed value.
			new_ni[get_index(row, col, M)] = numerator / denominator;
			
			// Update max_diff if applicable
			if (fabs(get(new_ni, row, col, M, N) - get(ni, row, col, M, N))) >
				*max_diff) {
				*max_diff = fabs(get(new_ni, row, col, M, N) - 
					get(ni, row, col, M, N));
			}
		}
	}
}

static double* update_p(Simulation* sim, double* max_diff) {
	double* new_p[sim -> M * sim -> N];
	
	for (int row = 0; row < sim -> M; row++) {
		for (int col = 1; col < (sim -> N) - 1; col++) {
			// S, T, U are intermediate variables to make computations easier
			double S = 0;
			double T = 0;
			
			S += E_R * get_left(sim -> p, row, col, sim -> M, sim -> N) / 
				(sim -> H * sim -> H);	// Left
			S += E_R * get_right(sim -> p, row, col, sim -> M, sim -> N) / 
				(sim -> H * sim -> H);	// Right
			S += E_R * get_top(sim -> p, row, col, sim -> M, sim -> N) / 
				(sim -> H * sim -> H);	// Top
			S += E_R * get_bottom(sim -> p, row, col, sim -> M, sim -> N) / 
				(sim -> H * sim -> H);	// Bottom
				
			T += E * get(sim -> ni_na, row, col, sim -> M, sim -> N) / E_0;
			T += E * get(sim -> ni_k, row, col, sim -> M, sim -> N) / E_0;
			T -= E * get(sim -> ni_cl, row, col, sim -> M, sim -> N) / E_0;
			T -= E * get(sim -> ni_a, row, col, sim -> M, sim -> N) / E_0;
			
			// Set element in new_p to computed value.
			new_p[get_index(row, col, sim -> M)] = S + T / (4 * E_R / 
				(sim -> H * sim -> H));
			
			// Update max_diff if applicable
			if (fabs(get(new_p, row, col, sim -> M, sim -> N) - 
				get(p, row, col, sim -> M, sim -> N))) > *max_diff) {
				*max_diff = fabs(get(new_p, row, col, sim -> M, sim -> N) - 
					get(p, row, col, sim -> M, sim -> N));
			}
		}
	} 
	
	return new_p;
}