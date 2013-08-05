function I = current(W, L, VGS_TOP, VSD)
% CURRENT Compute current through sensor of length L by integration
% Michelle Shu | July 24, 2013

Q = 1.60217657e-19; % Elementary charge in C
N_0 = 1.57e17;      % Charge carrier density (m^-2)
MOB = 0.1;          % m^2/s
VSAT = 6.5e4;       % m/s; lower bound in Meric paper
% VSAT = 5.5e5;     % m/s; upper bound in Meric paper

K = 1e-3;           % Grid size along V axis
V_values = (0.033 : K : VSD - 0.033); % Values of V up to Vsd for integration

% Include V0 adjustment
VGS_TOP_0 = 0;   % V
VGS_BACK_0 = 0;   % V
VGS_BACK = 0;       % V
C_BACK = 1.2e-4;    % F / m^2
V_0 = VGS_TOP_0 + (C_BACK ./ C_top_5nmEFF(V_values)) * (VGS_BACK_0 - VGS_BACK);
% V_0 = 0;

N = sqrt(N_0 ^ 2 + (C_top_5nmEFF(V_values) .* (VGS_TOP - V_values - V_0) / Q) .^ 2);
denominator = 1 + (MOB * VSD / (L * VSAT));

I = (Q * W * MOB / L) * trapz(V_values, N) / denominator;

end