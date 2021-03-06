function F = V_functionsB(V, I, W, H, VGS_TOP, VSD)
% V_FUNCTIONS Helper for GVOLTAGE - Generate matrix of all functions in our
% system.

    Q = 1.60217657e-19; % Elementary charge
    N_0 = 1.57e17;      % Charge carrier density (1/m^2)
    MOB = 0.1;          % m^2/s
    VSAT = 6.5e4;       % m/s; lower bound in Meric paper
    % VSAT = 5.5e5;     % m/s; upper bound in Meric paper

    % Constants for quantum capacitance modification
    % V_F = 1e6;                  % Fermi velocity (m/s)
    % DIRAC = 1.054571726e-34;    % Dirac constant (Js)

    F = zeros(numel(V) + 1, 1);
    F(numel(V)) = V(numel(V)) - VSD;
    F(numel(V) + 1) = V(1);

    for i = 1 : numel(V) - 1
        
        % Include V0 adjustment
        %VGS_TOP_0 = -0.5;   % V
        %VGS_BACK_0 = 150;   % V
        %VGS_BACK = 0;       % V
        %C_BACK = 1.2e-4;    % F / m^2
        %V_0 = VGS_TOP_0 + (C_BACK / C_top_1nmEFF(V(i))) * ...
        %     (VGS_BACK_0 - VGS_BACK);
        V_0 = 0;
        
        % Calculate charge carrier density based on Cedl
        N = sqrt(N_0^2 + (C_top_3nmEFF(V(i)) * (VGS_TOP - V(i) - V_0) / Q) ^2);
        % C_Q = 2 * (Q^2) * sqrt(N) / (DIRAC * V_F * sqrt(pi)); % quantum cap
        
        % Then use series combination of Cedl and Cq for remaining
        % calculations.
        % C = 1 / ((1 / C_top(V(i))) + (1 / C_Q));  
        % N = sqrt(N_0^2 + (C * (VGS_TOP - V(i) - V_0) / Q) ^2);
        
        if (i == 1)
            E = (V(2) / H); % hold V(1) at 0
        elseif (i == numel(V) - 1)
            E = (VSD - V(i)) / H;
        else
            E = (V(i + 1) - V(i)) / H; 
        end
        
        % New version that includes modification with VSAT;
        % Don't know if this will work...
        F(i) = I + (Q * W * N) * (MOB * E) / (1 + (MOB/VSAT * abs(E)));
        
        % The original version without modification with VSAT
        % F(i) = I + (Q * W * MOB * N * E);
    end
end