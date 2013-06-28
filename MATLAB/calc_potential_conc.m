% MPB_POTENTIAL_CONC_PARAMS Script
% Display how electric potential profile and ion concentration distribution
% varies with the following parameters:
% bulk concentration (0.01 - 0.3 M), relative permittivity (30 - 120),
% effective ion size (0.1 - 5 nm), valence (+/- 1 or +/- 2)

% Michelle Shu | June 25, 2013

% Constants
K = 1.3806488e-23;      % Boltzmann constant (J/K)
T = 293;                % Temperature (K)
E = 1.60217657e-19;     % Elementary charge (C)
N_A = 6.0221413e23;     % Avogadro's number (1/mol)

% Defaults
P_0 = 0.025;            % Applied potential
C_0 = 150;
E_R = 78.3;
EFF = 1e-9;
Z = 1;

% Options for parameter variance
C_0_options = [10 50 100 150 300]; % in mol/m^3 = M * 1e3
E_R_options = [30 50 78.3 100 120];
EFF_options = [1e-10 5e-10 1e-9 3e-9 5e-9];
Z_options = [1 2];

% -------------------------------------------------------------------------
% A. Bulk Concentration

pb_P_C_0 = zeros(numel(C_0_options), 301);
pb_Cpos_C_0 = zeros(numel(C_0_options), 301);
pb_Cneg_C_0 = zeros(numel(C_0_options), 301);
mpb_P_C_0 = zeros(numel(C_0_options), 301);
mpb_Cpos_C_0 = zeros(numel(C_0_options), 301);
mpb_Cneg_C_0 = zeros(numel(C_0_options), 301);

for i = 1 : numel(C_0_options)
    [pb_X, pb_P_C_0(i, :), ~] = potential_1d(P_0, C_0_options(i), E_R, Z);
    pb_Cpos_C_0(i, :) = (C_0_options(i) .* (exp((- E / (K * T)) .* ...
        pb_P_C_0(i, :))));
    pb_Cneg_C_0(i, :) = (C_0_options(i) .* (exp((E / (K * T)) .* ...
        pb_P_C_0(i, :))));
    
    [mpb_X, mpb_P_C_0(i, :), ~] = mpb_potential(P_0, C_0_options(i), ...
        E_R, EFF, Z);
    mpb_Cpos_C_0(i, :) = (C_0_options(i) .* (exp((- E / (K * T)) .* ...
        mpb_P_C_0(i, :))));
    mpb_Cneg_C_0(i, :) = (C_0_options(i) .* (exp((E / (K * T)) .* ...
        mpb_P_C_0(i, :))));
end



% ------------------------------------------------------------------------
% B. Permittivity

pb_P_E_R = zeros(numel(E_R_options), 301);
pb_Cpos_E_R = zeros(numel(E_R_options), 301);
pb_Cneg_E_R = zeros(numel(E_R_options), 301);
mpb_P_E_R = zeros(numel(E_R_options), 301);
mpb_Cpos_E_R = zeros(numel(E_R_options), 301);
mpb_Cneg_E_R = zeros(numel(E_R_options), 301);

for i = 1 : numel(E_R_options)
    [pb_X, pb_P_E_R(i, :), ~] = pb_potential(P_0, C_0, E_R_options(i), Z);
    pb_Cpos_E_R(i, :) = (C_0 .* (exp((- E / (K * T)) .* pb_P_E_R(i, :))));
    pb_Cneg_E_R(i, :) = (C_0 .* (exp((E / (K * T)) .* pb_P_E_R(i, :))));
    
    [mpb_X, mpb_P_E_R(i, :), ~] = mpb_potential(P_0, C_0, ...
        E_R_options(i), EFF, Z);
    mpb_Cpos_E_R(i, :) = (C_0 .* (exp((- E / (K * T)) .* mpb_P_E_R(i, :))));
    mpb_Cneg_E_R(i, :) = (C_0 .* (exp((E / (K * T)) .* mpb_P_E_R(i, :))));
end




%-------------------------------------------------------------------------
% C. Effective Ion Size

[pb_X, pb_P_EFF, ~] = pb_potential(P_0, C_0, E_R, Z);
pb_Cpos_EFF = (C_0 .* (exp((- E / (K * T)) .* pb_P_EFF)));
pb_Cneg_EFF = (C_0 .* (exp((E / (K * T)) .* pb_P_EFF)));


mpb_P_EFF = zeros(numel(EFF_options), 301);
mpb_Cpos_EFF = zeros(numel(EFF_options), 301);
mpb_Cneg_EFF = zeros(numel(EFF_options), 301);

for i = 1 : numel(EFF_options)
    [mpb_X, mpb_P_EFF(i, :), ~] = mpb_potential(P_0, C_0, ...
        E_R, EFF_options(i), Z);
    mpb_Cpos_EFF(i, :) = (C_0 .* (exp((- E / (K * T)) .* mpb_P_EFF(i, :))));
    mpb_Cneg_EFF(i, :) = (C_0 .* (exp((E / (K * T)) .* mpb_P_EFF(i, :))));
end




%--------------------------------------------------------------------------
% D. Valence

[pb_X, pb_P_Z1, ~] = pb_potential(P_0, C_0, E_R, Z_options(1));
[pb_X, pb_P_Z2, ~] = pb_potential(P_0, C_0, E_R, Z_options(2));
[mpb_X, mpb_P_Z1, ~] = mpb_potential(P_0, C_0, E_R, EFF, Z_options(1));
[mpb_X, mpb_P_Z2, ~] = mpb_potential(P_0, C_0, E_R, EFF, Z_options(2));

pb_Cpos_Z1 = (C_0 .* (exp((- E / (K * T)) .* pb_P_Z1)));
pb_Cneg_Z1 = (C_0 .* (exp((E / (K * T)) .* pb_P_Z1)));
pb_Cpos_Z2 = (C_0 .* (exp((- E / (K * T)) .* pb_P_Z2)));
pb_Cneg_Z2 = (C_0 .* (exp((E / (K * T)) .* pb_P_Z2)));

mpb_Cpos_Z1 = (C_0 .* (exp((- E / (K * T)) .* mpb_P_Z1)));
mpb_Cneg_Z1 = (C_0 .* (exp((E / (K * T)) .* mpb_P_Z1)));
mpb_Cpos_Z2 = (C_0 .* (exp((- E / (K * T)) .* mpb_P_Z2)));
mpb_Cneg_Z2 = (C_0 .* (exp((E / (K * T)) .* mpb_P_Z2)));


