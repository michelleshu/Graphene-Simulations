function [ X, P, R ] = potential_1d (P_0, Zi, Ci, E_R, EFF, MPB)
% PB_3P_DOUBLE_LAYER 1-D Simulation of electrical double layer using 
% Poisson-Boltzmann distribution and 3-point formula. Assumes point charge, 
% monovalent ions.

% USAGE 
%	Input parameters:
%	- P_0 : applied potential (V)
%	- Zi : ordered array of ion charge numbers (usually 1 or 2)
%	- Ci : ordered array of ion concentrations (in mol/m^3 -> M * 1e3)
%	- E_R : relative permittivity (78.3 for H2O)
%	- EFF : effective ion size (m, should overestimate at ~ 3e-8 m)
%   - MPB : 1 if use modified PB, 0 if use regular PB
%	Outputs:
%	- X : distance axis orthogonal to interface (m)
%	- P : matrix of potentials of points at x distance away
%	- R : matrix of second derivatives of potentials


%  PROCEDURE
%   1. Initial guess of potential distribution (P) is a decreasing linear 
%       function.
%   2. Use PB equation to compute second derivative of potential (R = P'')
%       with boundary conditions at current iteration by Poisson-Boltzmann
%   3. Derive new estimation for P (P_calc) from substituting R in
%       three point formula.
%   4. Check for convergence (P and P_calc differ by less than C). If not
%       converged, repeat all.
%   5. Adjust P closer to P_calc by ratio A.

%   Michelle Shu | June 17, 2013

% Constants
N_A = 6.0221413e23;     % Avogadro's number (1/mol)
E = 1.60217657e-19;     % Elementary charge (C)
E_0 = 8.854187817e-12;  % Vacuum permittivity (F/m)   
K = 1.3806488e-23;      % Boltzmann constant (J/K)
T = 293;                % Temperature (K)

H = 1e-10;              % Distance step size, determines resolution (m)
A = 1e-7 / (P_0 ^ 2);   % Initial potential function adjustment step ratio
                        %   (controls speed of convergence)
C = 1e-5;               % Convergence criterion: Max acceptable ratio of
                        %   (P_calc - P)/ P for an individual point
L = 3e-8;               % Limit for P to approach 0
G = 1e6;                % Steepness of initialization curve

% Initialization
X = (0 : H : L);
P = (P_0 / (1 - exp(-G * L))) * (exp(-G * X) - exp(-G * L));   % exp. guess
R = zeros(size(P));
P_calc = zeros(size(P));
P_calc(1) = P_0;   
done = false;

iter = 1;
while ~done 
    % (Modified) Poisson-Boltzmann equation
    for i = 1 : numel(P)
        R(i) = 0;
        for ion = 1 : numel(Zi)
            Z = Zi(ion);
            C_0 = Ci(ion);
            V = 2 * (EFF ^ 3) * C_0 * N_A;
            
            sinh_term1 = sinh(Z * E * P(i) / (K * T));
            R(i) = R(i) + (2 * Z * E * N_A * C_0 / (E_0 * E_R)) * ...
                sinh_term1;
            
            if (MPB == 1)
                sinh_term2 = sinh(Z * E * P(i) / (2 * K * T));
                R(i) = R(i) / (1 + (2 * V * (sinh_term2 ^ 2)));
            end
        end
    end
    
    % Three point formula
    for i = 2 : numel(P) - 1
        P_calc(i) = (P(i - 1) + P(i + 1) - R(i) * H^2) / 2;
    end
    max_difference = max(abs(P_calc - P) ./ P);
    % disp([num2str(max_difference) '   ' num2str(P_0)]);

    % Test for convergence
    if (max_difference) < C
        done = true;
    end
    
    P = P + (A * (P_calc - P));
    
    % Adjust step ratio to larger values gradually until it hits 0.01.
    % A starts at small value to avoid overshooting in cases with large P_0
    if (mod(iter, 100) == 0) && (A < 0.01)
        A = A * 10;
    end
    
    iter = iter + 1;
end

end

