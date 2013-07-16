/* SimulationGrid.h
 * Manage grid of ion density, potential, diffusion coefficient values
 * related to simulation.
 * Michelle Shu | July 16, 2013
 */

# include <stdlib.h>
# include <assert.h>
# include "Matrix.h"

/* ----- STRUCTURE DEFINITION ---------------------------------------------- */

typedef struct simulation Simulation;
struct simulation {
	int M;				// Number of rows
	int N;				// Number of columns
	int H;				// Grid size
	double* ni_na;		// Ion density
	double* ni_k;
	double* ni_cl;
	double* ni_a;
	double* d_na;		// Diffusion coefficient
	double* d_k;
	double* d_cl;
	double* d_a;
	double* p;			// Electrical potential
	int neuron_left;	// Leftmost column of neuron
	int neuron_right; 	// Rightmost column of neuron
	int neuron_top;		// Top row of neuron 
	int neuron_bottom;	// Bottom row of neuron
	int mem_squares;	// # of grid squares membrane takes up	
};

/* ----- PUBLIC FUNCTION PROTOTYPES ---------------------------------------- */

Simulation* simulation_new();
void simulation_destroy(Simulation* sim);