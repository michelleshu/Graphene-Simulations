V = (0 : 0.001 : 1);
N = zeros(size(V));
C_Q = zeros(size(V));
VGS_TOP = 0.4;

Q = 1.60217657e-19; % Elementary charge
N_0 = 1.57e17;      % Charge carrier density (m^-2)
MOB = 0.1;          % m^2/s
VSAT = 5.5e9;       % m/s; upper bound in Meric paper

V_F = 1e6;          % Fermi velocity (m/s)
DIRAC = 1.054571726e-34;    % Dirac constant (Js)


for i = 1 : numel(V) - 1
    % Calculate charge carrier density based on Cedl
    N(i) = sqrt(N_0^2 + (C_top(V(i)) * (VGS_TOP - V(i)) / Q) ^2);
    C_Q(i) = 2 * (Q^2) * sqrt(N(i)) / (DIRAC * V_F * sqrt(pi)); % quantum cap
end