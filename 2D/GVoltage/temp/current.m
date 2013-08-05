function I = current(W, L, VGS_TOP)
%CURRENT Compute current through sensor of length L by integration

Q = 1.60217657e-19; % Elementary charge
N_0 = 5e15; % Minimum carrier density
MOB = 0.055; % m^2 / s

K = 1e-3; % Grid size along V
V_values = (0 : K : VGS_TOP);

F1 = (C_top(V_values) .* (VGS_TOP - V_values) / Q) .^ 2;
F2 = sqrt((N_0^2) + F1);
I = Q * W * MOB/ L * trapz(F2, V_values); % Derive I by integration
end

