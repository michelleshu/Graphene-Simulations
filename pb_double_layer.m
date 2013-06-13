function pb_double_layer
% PB_DOUBLE_LAYER Simulation of electrical double layer potential due to
% Poisson-Boltzmann distribution. Assumes point charges, monovalent ions.

%   ** Set physical model parameters in PB_POTENTIAL.m **

%   Procedural details:
%   1. Initial guess of potential distribution (P) is a decreasing linear 
%       function.
%   2. Apply bvp4c and ode45 to PB_POTENTIAL to compute second derivative 
%       of potential (R = P'') with boundary conditions at current 
%       iteration by Poisson-Boltzmann equation.
%   3. Derive new estimation for P (P_CALC) from substituting R in
%       three point formula.
%   4. Check for convergence by measuring difference between P, P_CALC
%   

%   Michelle Shu | June 13, 2013

% Constants

H = 1e-10;              % Distance step size, determines resolution (m)
A = 0.1;                % Potential function adjustment step ratio
                        %   (controls speed of convergence)
C = 0.01;               % Convergence criterion: System reaches equilibrium
                        %   when mean percent error between P, P_CALC < C



end

