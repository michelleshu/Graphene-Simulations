% CHARGE_CAPACITANCE.M Script to compute and graph total charge and
% capacitance for Poisson Boltzmann and Modified Poisson Boltzmann models

% Across each row is points along x-axis for one potential function
% Down the rows are different functions for varying P0s.

% We sample 2e-8 / 1e-10 = 201 values for each applied potential P0.
% We try values from 10mV to 1V in 10mV increments for 100 values of P0.
% P will be a 100 x 201 matrix.
% R will be the 100 x 201 matrix holding the respective second derivatives
% of P.

% This will take a while so just save matrix P and R and do simulation in
% sections.

% Michelle Shu | June 21, 2013

P0_values = (0.01 : 0.01 : 1.00);  % Applied potential in volts
C_0 = 100;                         % Bulk concentration (mol/m^3)
E_R = 78.3;                        % Relative permittivity 
E_0 = 8.854187817e-12;

for index = 1 : 100
     disp(index);
     [pb_X, pb_P(index, :), pb_R(index, :)] = pb_potential(P0_values(index), C_0, E_R);
     [mpb_X, mpb_P(index, :), mpb_R(index, :)] = mpb_potential(P0_values(index), C_0, E_R);
end

% A. Total Charge

pb_Qt = charge_poisson(pb_X, pb_R, P0_values);
mpb_Qt = charge_poisson(mpb_X, mpb_R, P0_values);

semilogy(P0_values, mpb_Qt ./ 1e4, '-b', P0_values, pb_Qt ./ 1e4, '-r');
title(['Total Charge per Area v. Applied Voltage (P0)',...
    ' (z = 1, C0 = 0.1 M, er = 80, ion size = 0.3 nm)'], 'FontSize', 16);
xlabel('P0 (V)', 'FontSize', 16);
ylabel('Charge per Area (C / cm^2)', 'FontSize', 16);
legend('Modified Poisson Boltzmann', 'Poisson Boltzmann');

pause;

% B. Capacitance

pb_C = pb_Qt ./ P0_values;
mpb_C = mpb_Qt ./ P0_values;

semilogy(P0_values, mpb_C ./ 1e4, '-b', P0_values, pb_C ./ 1e4, '-r');
title(['Capacitance per Area v. Applied Potential (P0)',...
    '(z = 1, C0 = 0.1 M, er = 80, ion size = 0.3 nm)'], 'FontSize', 16);
xlabel('P0 (V)', 'FontSize', 16);
ylabel('Capacitance per Area (F / cm^2)', 'FontSize', 16);
legend('Modified Poisson Boltzmann', 'Poisson Boltzmann');

pause;

% C. Differential Capacitance

dV = P0_values(2) - P0_values(1);
pb_Cdiff = zeros(size(pb_Qt));
mpb_Cdiff = zeros(size(mpb_Qt));

for i = 1 : numel(pb_Qt) - 1
    pb_Cdiff(i + 1) = (pb_Qt(i + 1) - pb_Qt(i)) / dV;
    mpb_Cdiff(i + 1) = (mpb_Qt(i + 1) - mpb_Qt(i)) / dV;
end
semilogy(P0_values(1, 2 : end), mpb_Cdiff(1, 2 : end) ./ 1e4, '-b', ...
    P0_values(1, 2 : end), pb_Cdiff(1, 2 : end) ./ 1e4, '-r');
title(['Differential Capacitance per Area v. Applied Potential (P0)',...
    '(z = 1, C0 = 0.1 M, er = 80, ion size = 0.3 nm)'], 'FontSize', 16);
xlabel('P0 (V)', 'FontSize', 16);
ylabel('Differential Capacitance (F / cm ^2)', 'FontSize', 16);
legend('Modified Poisson Boltzmann', 'Poisson Boltzmann');
