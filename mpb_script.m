% Across each row is points along x-axis for one potential function
% Down the rows are different functions for varying P0s.

% We sample 2e-8 / 1e-10 = 201 values for each applied potential P0.
% We try values from 10mV to 1V in 10mV increments for 100 values of P0.
% P will be a 100 x 201 matrix.
% R will be the 100 x 201 matrix holding the respective second derivatives
% of P.

% This will take a while so just save matrix P and R and do simulation in
% sections.

P0_values = (0.01 : 0.01 : 1.00);    % Applied potential in volts
C_0 = 100;                         % Bulk concentration (mol/m^3)
E_R = 78.3;                        % Relative permittivity 
E_0 = 8.854187817e-12;

for index = 13 : 100
     disp(index);
     [pb_X, pb_P(index, :), pb_R(index, :)] = pb_potential(P0_values(index), C_0, E_R);
     [mpb_X, mpb_P(index, :), mpb_R(index, :)] = mpb_potential(P0_values(index), C_0, E_R);
end

% pb_Qt = charge_poisson(pb_X, pb_R, P0_values);
% mpb_Qt = charge_poisson(mpb_X, mpb_R, P0_values);

% semilogy(P0_values, mpb_Qt, '-b', P0_values, pb_Qt, '-g');
% title(['Total Charge (Q) v. Applied Voltage (P0)',...
%     '(z = 1, C0 = 0.1 M, er = 80)'], 'FontSize', 16);
% xlabel('P0 (V)', 'FontSize', 16);
% ylabel('Q (C)', 'FontSize', 16);