function [Qt] = charge_poisson (X, R, P0_values)
% PB_3P_CHARGE_POISSON Compute total charge contained in double layer by 
% integration of charge distribution derived from Poisson equation using 
% second derivative of potential. Graph relationship between total charge 
% and applied potential.

% Michelle Shu | June 18, 2013

E_0 = 8.854187817e-12;  % Vacuum permittivity (F/m)  
E_R = 79.3;             % Relative permittivity (80 for H2O)

Qt = zeros(size(P0_values));

for i = 1 : numel(P0_values)
    % Approximate integral by trapezoidal estimation for total charge Qt
    Qt(i) = trapz(X, R(i, :) .* (E_R * E_0));
end

plot(P0_values, Qt);
title(['Total Charge (Q) v. Total Voltage (P0)',...
    '(z = 1, C0 = 0.1 M, er = 80)'], 'FontSize', 16);
xlabel('P0 (V)', 'FontSize', 16);
ylabel('Q (C)', 'FontSize', 16);
end