% PB_3P_ION_CONCENTRATION Compute ion concentration profile near electrical
% double layer according to Boltzmann distribution

% Michelle Shu | June 18, 2013

% Constants
K = 1.3806488e-23;      % Boltzmann constant (J/K)
T = 293;                % Temperature (K)
E = 1.60217657e-19;     % Elementary charge (C)
N_A = 6.0221413e23;     % Avogadro's number (1/mol)

% Parameters
E_R = 80;
P_0 = 0.025;            % Surface potential [Range: 25 - 200 mV] (V)
C_0 = 100;              % Bulk concentration (mol/m^3) -> M * 1e3

[ X, P, ~ ] = pb_3p_double_layer(P_0, C_0, E_R);

Cpos = (C_0 .* (exp((- E / (K * T)) .* P))); % + ion concentration
Cneg = (C_0 .* (exp((E / (K * T)) .* P)));   % - ion concentration

semilogy(X, (Cpos ./ 1000), X, (Cneg ./ 1000));

title(['Ion Concentration v. Distance from Interface',...
    '(C0 = 0.1 M, P0 = 25 mV, er = 80)'], 'FontSize', 16);
xlabel('Distance (m)', 'FontSize', 16);
ylabel('Concentration (M)', 'FontSize', 16);
legend('z = +1', 'z = -1');