function [ X, P, R ] = pb_3p_double_layer (P_0, C_0, E_R)
% PB_3P_DOUBLE_LAYER Simulation of electrical double layer using 
% Poisson-Boltzmann distribution and 3-point formula. Assumes point charge, 
% monovalent ions.

%   Let P be the potential function (PSI) with respect to x (distance).
%   Let R be the second derivative of P (P'').

%   Procedure
%   1. Initial guess of potential distribution (P) is a decreasing linear 
%       function.
%   2. Use PB equation to compute second derivative of potential (R = P'')
%       with boundary conditions at current iteration by Poisson-Boltzmann
%   3. Derive new estimation for P (P_calc) from substituting R in
%       three point formula.
%   4. Check for convergence (P and P_calc differ by less than C). If not
%       converged, repeat all.
%   5. Adjust P closer to P_calc by ratio A.

%   Michelle Shu | Last modified on June 17, 2013

% Constants
N_A = 6.0221413e23;     % Avogadro's number (1/mol)
E = 1.60217657e-19;     % Elementary charge (C)
Z = 1;                  % Ion charge number
E_0 = 8.854187817e-12;  % Vacuum permittivity (F/m)   
K = 1.3806488e-23;      % Boltzmann constant (J/K)
T = 293;                % Temperature (K)

% Adjustable Parameters:
% Can set these by passing in as parameters
% P_0 = 0.025;            % Surface potential [Range: 25 - 200 mV] (V)
% C_0 = 100;              % Bulk concentration (mol/m^3) -> M * 1e3
% E_R = 78.3;               % Relative permittivity (80 for H2O)

H = 1e-10;              % Distance step size, determines resolution (m)
A = 0.1;                % Potential function adjustment step ratio
                        %   (controls speed of convergence)
C = 1e-6;               % Convergence criterion: Max acceptable ratio of
                        %   (P_calc - P)/ P for an individual point
L = 1e-8;               % Limit for P to approach 0
G = 1e6;                % Steepness of initialization curve

% Step 1. Initialization
X = (0 : H : L);
%P = (- P_0 / L) .* X + P_0;     % linear guess
P = (P_0 / (1 - exp(-G * L))) * (exp(-G * X) - exp(-G * L));   % exp. guess
R = zeros(size(P));
P_calc = zeros(size(P));
P_calc(1) = P_0;        % fix at P_0
done = false;

while ~done 
    % Step 2. Poisson-Boltzmann equation
    for i = 1 : numel(P)
        R(i) = (2 * Z * E * N_A * C_0 * sinh(Z * E * P(i) / (K * T))) / ...
            (E_R * E_0);
    end
    
    % Step 3. Three point formula
    for i = 2 : numel(P) - 1
        P_calc(i) = (P(i - 1) + P(i + 1) - R(i) * H^2) / 2;
    end
    max_difference = max(abs(P_calc - P) ./ P);
    disp(max_difference);

    % Step 4. Test for convergence
    if (max_difference) < C
        done = true;
    end
    
    % Step 5.
    P = P + (A * (P_calc - P));
end

semilogy(X, P);

end

