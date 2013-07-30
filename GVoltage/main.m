I = -3.6282086e-04; % Current (A)

W = 5e-6;           % Channel width (m)
L = 1e-5;           % Channel length (m)
H = 5e-8;           % Resolution (grid size along channel length

V_0 = 0.1;          % Applied potential
VGS_TOP = 0.4;      % Volts

[V, fval, X] = gvoltage(I, W, L, H, V_0, VGS_TOP);