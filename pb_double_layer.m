function [ R, P_calc ] = pb_double_layer
% PB_DOUBLE_LAYER Simulation of electrical double layer potential due to
% Poisson-Boltzmann distribution. Assumes point charges, monovalent ions.

%   Procedural details:
%   Method A
%   Solve Poisson-Boltzmann equation analytically 
%
%   Method B
%   1. Initial guess of potential distribution (P) is a decreasing linear 
%       function.
%   2. Use PB_POTENTIAL to compute second derivative of potential (R = P'')
%       with boundary conditions at current iteration by Poisson-Boltzmann
%   3. Derive new estimation for P (P_CALC) from substituting R in
%       three point formula.
%   4. Check for convergence
%   

%   Michelle Shu | June 13, 2013

% Parameters: P_0 is also defined in function RESIDUAL
P_0 = 0.1;              % Surface potential [Range: 25 - 200 mV] (V)
H = 1e-10;              % Distance step size, determines resolution (m)
A = 0.1;                % Potential function adjustment step ratio
                        %   (controls speed of convergence)
C = 0.01;               % Convergence criterion: System reaches equilibrium
                        %   when mean percent difference between 
                        

D = 2e-7;               % Limit for potential to reach 0                        
X = (0 : H : D);

solinit = bvpinit(X, [0 1]);
sol = bvp4c(@pb_function, @pb_boundary, solinit);
X = linspace(0, D);
P = deval(sol, X);
P

plot(X, P);

end

% ----------------------------------------------------------------------

function R = pb_function( X_range, P_init )
% PB_POTENTIAL Poisson-Boltzmann estimation of R = P''
%   P_init holds P(0) and P'(0) = - P_0 / D
%   dP/dx = P' and dP'/dx = R(2) = P'' by PB equation

% Constants
N_A = 6.0221413e23;     % Avogadro's number (1/mol)
E = 1.60217657e-19;     % Elementary charge (C)
Z = 1;                  % Ion charge number
K = 8.6173324e-5;       % Boltzmann constant (eV/K)
T = 20;                 % Temperature (K)

% Parameters:
E_R = 80;               % Relative permittivity (80 for H2O) 
C_0 = 1e-5;             % Bulk concentration (mol/L)

R = [P_init(2); - (Z * E * N_A * C_0 * sinh((Z * E * P_init(1)) / ...
    (K * T)) / E_R) ];
end

% ----------------------------------------------------------------------

function residual = pb_boundary (PA, PB)
% PB_BOUNDARY Boundary condition function for potential

P_0 = 0.1; 
residual = [(PA(1) - P_0) PB(1)];

end
