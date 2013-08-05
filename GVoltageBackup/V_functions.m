function F = V_functions(V, I, W, H, V_0, VGS_TOP)
% V_FUNCTIONS Helper for GVOLTAGE - Generate matrix of all functions in our
% system.

    Q = 1.60217657e-19; % Elementary charge
    N_0 = 5e15;         % Charge carrier density
    MOB = 0.055;        % m^2/s
    VSAT = 6.3e8;       % m/s; lower estimate in Meric paper

    V_F = 1e6;          % Fermi velocity (m/s)
    DIRAC = 1.054571726e-34;    % Dirac constant (Js)
    C_Q = 0.105;        % experimental value from Xia
    % C_Q = 2 * Q * Q * sqrt(N_0) / (DIRAC * V_F * sqrt(pi)); % Quantum C

    F = zeros(numel(V) + 1, 1);
    F(numel(V)) = V(numel(V)) - V_0;
    F(numel(V) + 1) = V(1);

    for i = 1 : numel(V) - 1
        C = 1 / ((1 / C_top(V(i))) + (1 / C_Q)); 
        E = (V(i + 1) - V(i)) / H;                      % electric field
        v_drift = MOB * E / (1 + (MOB * E / VSAT));     % drift velocity
        F(i) = I + (Q * W) * sqrt(N_0^2 + (C * ...
          (VGS_TOP - V(i)) / Q) ^ 2) * v_drift;
    end
end