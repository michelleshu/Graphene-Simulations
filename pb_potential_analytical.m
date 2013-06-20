function [ sol ] = pb_double_layer_analytical
% PB_DOUBLE_LAYER_ANALYTICAL Simulation of electrical double layer 
% potential due to Poisson-Boltzmann distribution. Assumes point charges, 
% monovalent ions.

%   Procedural details:
%   Solve Poisson-Boltzmann equation analytically using bvp4c

%   Michelle Shu | June 13, 2013

% Parameters:
P_0 = 0.1;              % Surface potential [Range: 25 - 200 mV] (V)
H = 1e-10;              % Distance step size, determines resolution (m)
A = 0.1;                % Potential function adjustment step ratio
                        %   (controls speed of convergence)
C = 0.01;               % Convergence criterion: System reaches equilibrium
                        %   when mean percent difference between 
                        

D = 1e-5;               % Limit for potential to reach 0                   

solinit = bvpinit(linspace(0, D, 2), @pb_init);
sol = bvp4c(@pb_ode, @pb_bc, solinit);
plot(sol.x, sol.y);

end

% ----------------------------------------------------------------------

function dPdx = pb_ode( X, P )
% PB_POTENTIAL Poisson-Boltzmann estimation of R = P''

%   P vector holds P(0) and P'(0)
%   dP/dx = P' and dP'/dx = R(2) = P'' by PB equation

% Constants
N_A = 6.0221413e23;     % Avogadro's number (1/mol)
E = 1.60217657e-19;     % Elementary charge (C)
Z = 1;                  % Ion charge number
K = 1.3806488e-23;      % Boltzmann constant (J/K)
T = 20;                 % Temperature (K)

% Parameters:
E_R = 80;               % Relative permittivity (80 for H2O) 
C_0 = 1e-2;             % Bulk concentration (mol/m^3)

dPdx = [P(2); (2 * Z * E * N_A * C_0 * sinh((Z * E * P(1)) / ...
    (K * T)) / E_R) ];
end

% ----------------------------------------------------------------------

function residual = pb_bc (PA, PB)
% PB_BOUNDARY Boundary condition function for potential

P_0 = 0.1; 
residual = [(PA(1) - P_0); PB(1)];

end

% ----------------------------------------------------------------------

function v = pb_init (X)
% PB_INIT Initial guess function for potential
P_0 = 0.1;
D = 1e-5;
v = [((- P_0 / D) * X + P_0) ; (- P_0 / D)];

end
