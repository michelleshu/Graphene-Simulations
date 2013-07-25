/* SimulationGrid.c
 * Manage grid of ion density, potential, diffusion coefficient values
 * related to simulation.
 * Michelle Shu | July 16, 2013
 */

# include "SimulationGrid.h"

const double N_A  = 6.0221413e23;		// Avogadro's number (1/mol)

/* ----- MODEL PARAMETERS -------------------------------------------------- */

/*// TOY CASE
const float GRID_WIDTH		= 10;     
const float GRID_HEIGHT    	= 10;     
const float NEURON_WIDTH	= 6; 
const float NEURON_HEIGHT	= 6;
const float NEURON_DISP		= 2;		// Displacement from interface
const float MEM_THICKNESS	= 1;		// Thickness of cell membrane
const float H   			= 1;	    // Size of grid square (resolution)

// Applied Potential in V
const double P_0			= 0.025;

// Ion Concentrations in mol/m^3
const float IC_CONC_NA		= 18;		// Intracellular
const float IC_CONC_K		= 135;
const float IC_CONC_CL		= 7;
const float IC_CONC_A		= 74;		// A- is generic organic anion
const float EC_CONC_NA		= 145;		// Extracellular
const float EC_CONC_K		= 3;
const float EC_CONC_CL		= 20;
const float EC_CONC_A		= 13;

// Diffusion Constants in m^2/s
const float D_NA			= 1;
const float D_K				= 1;
const float D_CL			= 1;
const float D_A				= 1;  // (Glutamate)

const float REST_MEM_D_NA	= 1e-6;  // D in membrane during resting state
const float REST_MEM_D_K	= 1;
const float REST_MEM_D_CL	= 1;
const float REST_MEM_D_A	= 1;*/

// Dimensions of simulation grid and neuron (sensor height = grid height) in m
const float GRID_WIDTH		= 1.1e-6;     
const float GRID_HEIGHT    	= 1.1e-6;     
const float NEURON_WIDTH	= 1e-6; 
const float NEURON_HEIGHT	= 1e-6;
const float NEURON_DISP		= 4e-8;		// Displacement from interface
const float MEM_THICKNESS	= 1e-8;		// Thickness of cell membrane
const float H   			= 1e-8;	    // Size of grid square (resolution)

// Applied Potential in V
const double P_0			= 0.2;

// Ion Concentrations in mol/m^3
const float IC_CONC_NA		= 18;		// Intracellular
const float IC_CONC_K		= 135;
const float IC_CONC_CL		= 7;
const float IC_CONC_A		= 74;		// A- is generic organic anion
const float EC_CONC_NA		= 145;		// Extracellular
const float EC_CONC_K		= 3;
const float EC_CONC_CL		= 20;
const float EC_CONC_A		= 13;

// Diffusion Constants in m^2/s
const float D_NA			= 1.33e-9;
const float D_K				= 1.96e-9;
const float D_CL			= 2.03e-9;
const float D_A				= 0.76e-9;  // (Glutamate)

const float REST_MEM_D_NA	= 0.04e-9;  // D in membrane during resting state
const float REST_MEM_D_K	= 1.00e-9;
const float REST_MEM_D_CL	= 0.45e-9;
const float REST_MEM_D_A	= 0.10e-9;

/* ----- PRIVATE FUNCTION PROTOTYPES --------------------------------------- */

static Simulation* simulation_allocate();
static void simulation_initialize(Simulation* sim);
static void set_extracellular(Simulation* sim, int i);
static void set_membrane(Simulation* sim, int i);
static void set_intracellular(Simulation* sim, int i);

/* ----- PUBLIC FUNCTIONS -------------------------------------------------- */

Simulation* simulation_new() {
	Simulation* sim = simulation_allocate();
	simulation_initialize(sim);
	return sim;
}

void simulation_destroy(Simulation* sim) {	
	free(sim -> ni_na);
	free(sim -> ni_k);
	free(sim -> ni_cl);
	free(sim -> ni_a);
	free(sim -> d_na);
	free(sim -> p);	
	free(sim);
}

/* ----- PRIVATE FUNCTIONS ------------------------------------------------- */

// Allocate space for all variables in simulation
static Simulation* simulation_allocate() {	
	Simulation* sim = (Simulation *) malloc(sizeof(Simulation));
	
	// Calculate number of rows and columns
	int M = (int) (GRID_HEIGHT / H);
	int N = (int) (GRID_WIDTH / H);
	sim -> M = M;
	sim -> N = N;
	sim -> H = H;
	
	sim -> ni_na = malloc(M * N * sizeof(double));
	sim -> ni_k = malloc(M * N * sizeof(double));
	sim -> ni_cl = malloc(M * N * sizeof(double));
	sim -> ni_a = malloc(M * N * sizeof(double));
	sim -> d_na = malloc(M * N * sizeof(double));
	sim -> d_k = malloc(M * N * sizeof(double));
	sim -> d_cl = malloc(M * N * sizeof(double));
	sim -> d_a = malloc(M * N * sizeof(double));
	sim -> p = malloc(M * N * sizeof(double));
	
	assert(sim);
	assert(sim -> ni_na);
	assert(sim -> ni_k);
	assert(sim -> ni_cl);
	assert(sim -> ni_a);
	assert(sim -> d_na);
	assert(sim -> d_k);
	assert(sim -> d_cl);
	assert(sim -> d_a);
	assert(sim -> p);

	return sim;
}

static void simulation_initialize(Simulation* sim) {
	// Record neuron placement in the grid
	sim -> neuron_left = (int) (NEURON_DISP / H);
	sim -> neuron_right = (int) ((NEURON_DISP + NEURON_WIDTH) / H);
	sim -> neuron_top = (int) ((GRID_HEIGHT - NEURON_HEIGHT) / (2 * H));
	sim -> neuron_bottom = (int) (sim -> neuron_top + (NEURON_HEIGHT / H));
	sim -> mem_squares = (int) (MEM_THICKNESS / H);
	
	// Start by setting every grid entry to extracellular conditions
	for (int i = 0; i < (sim -> M * sim -> N); i++) {
		set_extracellular(sim, i);
	}
	
	// Set membrane conditions
	// Top border
	for (int row = sim -> neuron_top; row < (sim -> neuron_top + 
		sim -> mem_squares); row++) {
		for (int col = sim -> neuron_left; col < sim -> neuron_right; 
			col++) {
			set_membrane(sim, get_index(row, col, sim -> M));
		}
	}
	// Left border
	for (int row = sim -> neuron_top; row < sim -> neuron_bottom; row++) {
		for (int col = sim -> neuron_left; col < (sim -> neuron_left + 
			sim -> mem_squares); col++) {
			set_membrane(sim, get_index(row, col, sim -> M));
		}
	}
	// Bottom border
	for (int row = (sim -> neuron_bottom - sim -> mem_squares); 
		row < sim -> neuron_bottom; row++) {
		for (int col = sim -> neuron_left; col < sim -> neuron_right;
			col++) {
			set_membrane(sim, get_index(row, col, sim -> M));	
		}
	}
	// Right border
	for (int row = sim -> neuron_top; row < sim -> neuron_bottom; row++) {
		for (int col = (sim -> neuron_right - sim -> mem_squares); 
			col < sim -> neuron_right; col++) {
			set_membrane(sim, get_index(row, col, sim -> M));	
		}
	}
	
	// Set intracellular conditions
	for (int row = (sim -> neuron_top + sim -> mem_squares);
		row < (sim -> neuron_bottom - sim -> mem_squares); row++) {
		for (int col = (sim -> neuron_left + sim -> mem_squares);
			col < (sim -> neuron_right - sim -> mem_squares); col++) {
			set_intracellular(sim, get_index(row, col, sim -> M));	
		}	
	}
	
	// Set potential
	for (int col = 0; col < sim -> N; col++) {
		for (int row = 0; row < sim -> M; row++) {
			sim -> p[get_index(row, col, sim -> M)] = (- P_0 / (GRID_WIDTH - H))
				* (col * H) + P_0;
		}
	}
	
	// Set concentration boundary condition at interface
	for (int row = 0; row < sim -> M; row++) {
		sim -> ni_na[get_index(row, 0, sim -> M)] = 0.0;
		sim -> ni_k[get_index(row, 0, sim -> M)] = 0.0;
		sim -> ni_cl[get_index(row, 0, sim -> M)] = 0.0;
		sim -> ni_a[get_index(row, 0, sim -> M)] = 0.0;
	}
}

/* Set cell at index i in simulation grid to extracellular conditions */
static void set_extracellular(Simulation* sim, int i) {
	sim -> ni_na[i] = EC_CONC_NA * N_A * H * H; // convert to # density
	sim -> ni_k[i] = EC_CONC_K * N_A * H * H;
	sim -> ni_cl[i] = EC_CONC_CL * N_A * H * H;
	sim -> ni_a[i] = EC_CONC_A * N_A * H * H;
	
	sim -> d_na[i] = D_NA;
	sim -> d_k[i] = D_K;
	sim -> d_cl[i] = D_CL;
	sim -> d_a[i] = D_A;
}

/* Set cell at index i in simulation grid to membrane conditions */
static void set_membrane(Simulation* sim, int i) {
	sim -> ni_na[i] = IC_CONC_NA * N_A * H * H; // Assume same ion concentration as inside cell
	sim -> ni_k[i] = IC_CONC_K * N_A * H * H;
	sim -> ni_cl[i] = IC_CONC_CL * N_A * H * H;
	sim -> ni_a[i] = IC_CONC_A * N_A * H * H;
	
	sim -> d_na[i] = REST_MEM_D_NA;
	sim -> d_k[i] = REST_MEM_D_K;
	sim -> d_cl[i] = REST_MEM_D_CL;
	sim -> d_a[i] = REST_MEM_D_A;
}

static void set_intracellular(Simulation* sim, int i) {
	sim -> ni_na[i] = IC_CONC_NA * N_A * H * H; // convert to # density
	sim -> ni_k[i] = IC_CONC_K * N_A * H * H;
	sim -> ni_cl[i] = IC_CONC_CL * N_A * H * H;
	sim -> ni_a[i] = IC_CONC_A * N_A * H * H;
	
	sim -> d_na[i] = D_NA;
	sim -> d_k[i] = D_K;
	sim -> d_cl[i] = D_CL;
	sim -> d_a[i] = D_A;
}