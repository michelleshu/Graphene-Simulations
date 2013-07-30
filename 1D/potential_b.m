N_A = 6.0221413e23;     % Avogadro's number (1/mol)
E = 1.60217657e-19;     % Elementary charge (C)
E_0 = 8.854187817e-12;  % Vacuum permittivity (F/m)
E_R = 78.3;
K = 1.3806488e-23;      % Boltzmann constant (J/K)
T = 293;                % Temperature (K)

C = 100;                % Concentration 0.1 mol / m^3
N = C * N_A;            % atoms / m^3
V0 = 1;              % 50 mV

X = (0 : 1e-10 : 2e-8);

KAPPA = sqrt(2 * N * (E^2) / (E_0 * E_R * K * T));

inner = tanh(E * V0 / (4 * K * T)) * exp(-KAPPA .* X);
potential = (4 * K * T / E) .* atanh(inner);