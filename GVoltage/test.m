Q = 1.60217657e-19; % Elementary charge
I = -2.418576e-4;   % Current (A)

W = 5e-6;           % Channel width
H = 1e-8;           % Grid size along sensor length

VGS_TOP = 0.4;      % Volts
N_0 = 5e15;         % Minimum carrier density
MOB = 0.055;        % m^2 / s

residual = zeros(size(V));

for i = 1 : numel(V) - 1
    residual(i) = I + (Q * W * MOB) * sqrt(N_0^2 + (C_top(V(i)) * ...
      (VGS_TOP - V(i)) / Q) ^ 2) * (V(i + 1) - V(i)) / (H);
end