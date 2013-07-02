function [X, Y, P, R] = potential_2d (P_0, Zi, Ci, E_R, EFF, MPB)
% POTENTIAL_2D   Evaluate potential distribution of 2-dimensional matrix

% USAGE 
%	Input parameters:
%	- P_0 : applied potential (V)
%	- Zi : ordered array of ion charge numbers (usually 1 or 2)
%	- Ci : ordered array of ion concentrations (in mol/m^3 -> M * 1e3)
%	- E_R : relative permittivity (78.3 for H2O)
%	- EFF : effective ion size (m, should overestimate at ~ 3e-8 m)
%   - MPB : 1 if use modified PB, 0 if use regular PB
%	Outputs:
%	- X : distance axis orthogonal to plate (m)
% 	- Y : distance axis parallel to plate (m)
%	- P : matrix of potentials at points (x, y)
%	- R : matrix of second derivatives of potentials at points (x, y)

%	Michelle Shu | June 23, 2013

% Constants
N_A = 6.0221413e23;     % Avogadro's number (1/mol)
E = 1.60217657e-19;     % Elementary charge (C)
E_0 = 8.854187817e-12;  % Vacuum permittivity (F/m)   
K = 1.3806488e-23;      % Boltzmann constant (J/K)
T = 293;                % Temperature (K)

% Adjustable Parameters
H = 1e-10;               % Distance step size, determines resolution (m)
A = 0.1;               % Initial potential function adjustment step ratio
                         %   (controls speed of convergence)
CONV = 1e-4;             % Convergence criterion: Max acceptable ratio of
                         %   (P_calc - P)/ P for an individual point
D = 2e-8;                % Limit for P to approach 0 (Max X value)
L = 1e-8;                % Length of interface
G = 1e6;                 % Steepness of initialization curve


% Initialization
X = (0 : H : D);
Y = (0 : H : L);
P = zeros(numel(Y), numel(X));
for i = 1 : numel(Y)
    P(i, :) = (P_0 / (1 - exp(-G * D))) * (exp(-G * X) - exp(-G * D));
end;
    
R = zeros(size(P));
P_calc = zeros(size(P));
P_calc(:, 1) = P_0;
done = false;

while ~done 
    % (Modified) Poisson-Boltzmann equation
    for row = 1 : numel(Y)
        for col = 1 : numel(X)
            R(row, col) = 0;
            for ion = 1 : numel(Zi)
                Z = Zi(ion);
                C_0 = Ci(ion);
                V = 2 * (EFF ^ 3) * C_0 * N_A; 
            
                sinh_term1 = sinh(Z * E * P(row, col) / (K * T));
                R(row, col) = R(row, col) + ...
                    (2 * Z * E * N_A * C_0/ (E_0 * E_R)) * sinh_term1;
                
                if (MPB == 1)
                    sinh_term2 = sinh(Z * E * P(row, col) / (2 * K * T));
                    R(row, col) = R(row, col) / ...
                    (1 + (2 * V * (sinh_term2 ^ 2)));
                end
            end
        end
    end
    
    % Update P based on neighbors
    for row = 1 : numel(Y)
        for col = 2 : numel(X) - 1
            P_calc(row, col) = ...
                (p_neighbors(P, col, row, numel(X), numel(Y)) - ...
                (R(row, col) * H * H)) / 4;
        end
    end
    
    max_difference = max(max(abs(P_calc - P) ./ P));
    %disp([num2str(max_difference) '   ' num2str(P_0)]);

    % Step 4. Test for convergence
    if (max_difference) < CONV
        done = true;
    end
    
    % Step 5.
    P = P + (A * (P_calc - P));

end
end