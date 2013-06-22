function [Qt] = charge_boltzmann (X, P, P0_values)
% PB_3P_CHARGE_BOLTZMANN Compute total charge contained in double layer by 
% first computing + and - ion concentration distribution in EDL with same 
% method as PB_3P_ION_CONCENTRATION, then converting to charge distribution
% and integrating to compute charge. Graph relationship between total 
% charge and applied potential.

% Michelle Shu | June 18, 2013

% Constants
K = 1.3806488e-23;      % Boltzmann constant (J/K)
T = 293;                % Temperature (K)
E = 1.60217657e-19;     % Elementary charge (C)
N_A = 6.0221413e23;     % Avogadro's number (1/mol)
EFF = 1e-10;               % Effective ion size (m) 

% Parameters
                        % Surface potential [Range: 25 - 200 mV] (V)
C_0 = 100;              % Bulk concentration (mol/m^3) -> M * 1e3

Qt = zeros(size(P0_values));
V = (2 * C_0 * N_A) * (EFF ^ 3);


for i = 1 : numel(P0_values)
    
    % Denominator for modification that includes steric effects
    steric_mod_denom = 1 + ((2 * V) .* ((sinh((E/(2 * K * T)) .* P)) .^ 2));
    
    Cpos = (C_0 .* (exp((- E / (K * T)) .* P))) ./ steric_mod_denom; 
                                                % + ion concentration
    Cneg = (C_0 .* (exp((E / (K * T)) .* P))) ./ steric_mod_denom;   
                                                % - ion concentration
    q = (Cneg .* (E * N_A)) - (Cpos .* (E * N_A)); % charge distribution
    
    % Approximate integral by trapezoidal estimation for total charge Qt
    Qt(i) = trapz(X, q);
end

plot(P0_values, Qt);
title(['Total Charge per Area v. Total Voltage (P0)',...
    ' (z = 1, C0 = 0.1 M, er = 80)'], 'FontSize', 16);
xlabel('P0 (V)', 'FontSize', 16);
ylabel('Charge per Area (C / m^2)', 'FontSize', 16);

end