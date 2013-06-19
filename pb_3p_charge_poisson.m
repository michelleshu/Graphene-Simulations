function Qt = pb_3p_charge_poisson
% PB_3P_CHARGE_POISSON Compute total charge contained in double layer by 
% integration of charge distribution derived from Poisson equation using 
% second derivative of potential. Graph relationship between total charge
% and applied potential.

% Michelle Shu | June 18, 2013

% Constants
E_0 = 8.854187817e-12;  % Vacuum permittivity (F/m)  

% Parameters
P_0_options = [0.025 : 0.005 : 0.200];            
                        % Surface potential [Range: 25 - 200 mV] (V)
C_0 = 100;              % Bulk concentration (mol/m^3) -> M * 1e3
E_R = 80;               % Relative permittivity (80 for H2O)

Qt = zeros(size(P_0_options));

for i = 1 : numel(P_0_options)
    [ X, ~, R ] = pb_3p_double_layer(P_0_options(i), C_0, E_R);
    q = R .* (E_R * E_0);
    % Approximate integral by trapezoidal estimation for total charge Qt
    Qt(i) = trapz(X, q);
end

semilogy(P_0_options, Qt);
title(['Total Charge (Q) v. Total Voltage (P0)',...
    '(z = 1, C0 = 0.1 M, er = 80)'], 'FontSize', 16);
xlabel('P0 (V)', 'FontSize', 16);
ylabel('Q (C)', 'FontSize', 16);
end