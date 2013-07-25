function F = V_functions(V)
% V_FUNCTIONS Helper for GVOLTAGE - Generate matrix of all functions in our
% system.

W = 5e-6;           % Channel width
Q = 1.60217657e-19; % Elementary charge
V_0 = 0.1;          % Applied potential
N_0 = 5e15;         % Minimum carrier density
MOB = 0.055;        % m^2 / s
H = 1e-7;           % Grid size along sensor length
I = -2.418576e-4;   % Current (A)
VGS_TOP = 0.4;      % Volts

F = zeros(numel(V), 1);
F(1) = V(1);
F(numel(V)) = V(numel(V)) - V_0;

for i = 2 : numel(V) - 1
    F(i) = I + (Q * W * MOB) * sqrt(N_0^2 + (C_top(V(i)) * ...
      (VGS_TOP - V(i)) / Q) ^ 2) * (V(i + 1) - V(i - 1)) / (2 * H);
end