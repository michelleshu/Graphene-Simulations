function [ X, P ] = pb_3p_double_layer
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

%   Michelle Shu | June 14, 2013

% Constants
N_A = 6.0221413e23;     % Avogadro's number (1/mol)
E = 1.60217657e-19;     % Elementary charge (C)
Z = 1;                  % Ion charge number
E_R = 80;               % Relative permittivity (80 for H2O) 
K = 1.3806488e-23;      % Boltzmann constant (J/K)
T = 20;                 % Temperature (K)

% Adjustable Parameters:
P_0 = 0.1;              % Surface potential [Range: 25 - 200 mV] (V)
C_0 = 1e-2;             % Bulk concentration (mol/m^3)
H = 1e-8;               % Distance step size, determines resolution (m)
A = 1e-11;              % Potential function adjustment step ratio
                        %   (controls speed of convergence)
C = 1;                  % Convergence criterion: Max acceptable ratio of
                        %   (P_calc - P)/ P
L = 1e-5;               % Limit for P to approach 0

% Step 1. (Initialization)
X = (0 : H : L);
P = (- P_0 / L) .* X + P_0;
R = zeros(size(P));
P_calc = zeros(size(P));
P_calc(1) = P_0;        % fix at P_0
done = false;

while ~done
    % Step 2.
    for i = 1 : numel(P)
        R(i) = (2 * Z * E * N_A * C_0 * sinh((Z * E * P(i)) / (K * T)) ...
            / E_R);     % Poisson-Boltzmann equation
    end
    
    % Step 3.
    residual_sum = 0;   % Compute these to prepare for step 4
    P_sum = 0;
    for i = 2 : numel(P) - 1
        P_calc(i) = (P(i - 1) + P(i + 1) - R(i) * H^2) / 2; % 3 point form.
        residual_sum = residual_sum + abs(P_calc(i) - P(i));
        P_sum = P_sum + P(i);
    end

    % Step 4.
    disp(residual_sum / P_sum);
    if (residual_sum / P_sum) < C
        done = true;
    end
    
    % Step 5.
    P = P + (A * (P_calc - P));
end

plot(X, P);

end

