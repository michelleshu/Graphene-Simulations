% CALC_CHARGE_CAP.M
% Compute functions for charge, capacitance and differential capacitance
% due to Modified Poisson Boltzmann model while varying the following 
% parameters: 
% bulk concentration (0.01 - 0.3 M), relative permittivity (30 - 120),
% effective ion size (0.1 - 5 nm), valence (+/- 1 or +/- 2)

% Michelle Shu | June 23, 2013

% Defaults
C_0 = 150;  % mol/m^3
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

pb_Qt_C_0 = zeros(numel(C_0_options), numel(P0_values));
mpb_Qt_C_0 = zeros(numel(C_0_options), numel(P0_values));

for i = 1 : numel(C_0_options)
    pb_P = zeros(numel(P0_values), 301);
    pb_R = zeros(numel(P0_values), 301);
    mpb_P = zeros(numel(P0_values), 301);
    mpb_R = zeros(numel(P0_values), 301);
    
    for j = 1 : numel(P0_values)
        disp([num2str(i) '   ' num2str(j)]);
        [pb_X, pb_P(j, :), pb_R(j, :)] = potential_1d(P0_values(j), ...
            Z, C_0_options(i), E_R, EFF, 0);
        [mpb_X, mpb_P(j, :), mpb_R(j, :)] = potential_1d(P0_values(j),...
            Z, C_0_options(i), E_R, EFF, 1);
    end
    
    pb_Qt_C_0(i, :) = charge_poisson(pb_X, pb_R, P0_values);
    mpb_Qt_C_0(i, :) = charge_poisson(mpb_X, mpb_R, P0_values);
end


pb_C_C_0 = zeros(size(pb_Qt_C_0));
mpb_C_C_0 = zeros(size(mpb_Qt_C_0));

for i = 1 : 5
    pb_C_C_0(i, :) = pb_Qt_C_0(i, :) ./ P0_values;
    mpb_C_C_0(i, :) = mpb_Qt_C_0(i, :) ./ P0_values;
end

pb_Cdiff_C_0 = zeros(size(pb_Qt_C_0));
mpb_Cdiff_C_0 = zeros(size(mpb_Qt_C_0));

for i = 1 : 5
    for j = 1 : numel(P0_values) - 1
        dV = P0_values(j + 1) - P0_values(j);
        pb_Cdiff_C_0(i, j + 1) = (pb_Qt_C_0(i, j + 1) - pb_Qt_C_0(i, j)) / dV;
        mpb_Cdiff_C_0(i, j + 1) = (mpb_Qt_C_0(i, j + 1) - mpb_Qt_C_0(i, j)) / dV;
    end
end


% -------------------------------------------------------------------------
% B. Relative Permittivity

pb_Qt_E_R = zeros(numel(E_R_options), numel(P0_values));
mpb_Qt_E_R = zeros(numel(E_R_options), numel(P0_values));

for i = 1 : numel(E_R_options)
    pb_P = zeros(numel(P0_values), 301);
    pb_R = zeros(numel(P0_values), 301);
    mpb_P = zeros(numel(P0_values), 301);
    mpb_R = zeros(numel(P0_values), 301);
    
    for j = 1 : numel(P0_values)
        disp([num2str(i) '   ' num2str(j)]);
        [pb_X, pb_P(j, :), pb_R(j, :)] = potential_1d(P0_values(j), ...
            Z, C_0, E_R_options(i), EFF, 0);
        [mpb_X, mpb_P(j, :), mpb_R(j, :)] = potential_1d(P0_values(j),...
            Z, C_0, E_R_options(i), EFF, 1);
    end
    
    pb_Qt_E_R(i, :) = charge_poisson(pb_X, pb_R, P0_values);
    mpb_Qt_E_R(i, :) = charge_poisson(mpb_X, mpb_R, P0_values);
end

pb_C_E_R = zeros(size(pb_Qt_C_0));
mpb_C_E_R = zeros(size(mpb_Qt_C_0));

for i = 1 : 5
    pb_C_E_R(i, :) = pb_Qt_E_R(i, :) ./ P0_values;
    mpb_C_E_R(i, :) = mpb_Qt_E_R(i, :) ./ P0_values;
end

pb_Cdiff_E_R = zeros(size(pb_Qt_E_R));
mpb_Cdiff_E_R = zeros(size(mpb_Qt_E_R));

for i = 1 : 5
    for j = 1 : numel(P0_values) - 1
        dV = P0_values(j + 1) - P0_values(j);
        pb_Cdiff_E_R(i, j + 1) = (pb_Qt_E_R(i, j + 1) - pb_Qt_E_R(i, j)) / dV;
        mpb_Cdiff_E_R(i, j + 1) = (mpb_Qt_E_R(i, j + 1) - mpb_Qt_E_R(i, j)) / dV;
    end
end



% ------------------------------------------------------------------------
% C. Effective Ion Size

mpb_Qt_EFF = zeros(numel(EFF_options), numel(P0_values));

for j = 1 : numel(P0_values)
    disp([num2str(i) '   ' num2str(j)]);
    [pb_X, pb_P(j, :), pb_R(j, :)] = potential_1d(P0_values(j), ...
        Z, C_0, E_R, EFF, 0);
end
pb_Qt_EFF = charge_poisson(pb_X, pb_R, P0_values);

for i = 1 : numel(EFF_options)
    mpb_P = zeros(numel(P0_values), 301);
    mpb_R = zeros(numel(P0_values), 301);
    
    for j = 1 : numel(P0_values)
        disp([num2str(i) '   ' num2str(j)]);
        [mpb_X, mpb_P(j, :), mpb_R(j, :)] = potential_1d(P0_values(j),...
            Z, C_0, E_R, EFF_options(i), 1);
    end
    
    mpb_Qt_EFF(i, :) = charge_poisson(mpb_X, mpb_R, P0_values);
end


pb_C_EFF = pb_Qt_EFF(1, :) ./ P0_values;

mpb_C_EFF = zeros(size(mpb_Qt_EFF));
for i = 1 : 5
    mpb_C_EFF(i, :) = mpb_Qt_EFF(i, :) ./ P0_values;
end

dV = P0_values(2) - P0_values(1);

pb_Cdiff_EFF = zeros(size(pb_Qt_EFF));
mpb_Cdiff_EFF = zeros(size(mpb_Qt_EFF));

for j = 1 : numel(P0_values) - 1
    pb_Cdiff_EFF(1, j + 1) = (pb_Qt_EFF(1, j + 1) - pb_Qt_EFF(1, j)) / dV;
end

for i = 1 : 5
    for j = 1 : numel(P0_values) - 1
        mpb_Cdiff_EFF(i, j + 1) = (mpb_Qt_EFF(i, j + 1) - mpb_Qt_EFF(i, j)) / dV;
    end
end

% -------------------------------------------------------------------------
% D. Valence

pb_Qt_Z = zeros(numel(Z_options), numel(P0_values));
mpb_Qt_Z = zeros(numel(Z_options), numel(P0_values));
 
for i = 1 : numel(Z_options)
    pb_P = zeros(numel(P0_values), 301);
    pb_R = zeros(numel(P0_values), 301);
    mpb_P = zeros(numel(P0_values), 301);
    mpb_R = zeros(numel(P0_values), 301);
    
    for j = 1 : numel(P0_values)
        disp([num2str(i) '   ' num2str(j)]);
        [pb_X, pb_P(j, :), pb_R(j, :)] = potential_1d(P0_values(j), ...
            Z_options(i), C_0, E_R, EFF, 0);
        [mpb_X, mpb_P(j, :), mpb_R(j, :)] = potential_1d(P0_values(j),...
            Z_options(i), C_0, E_R, EFF, 1);
    end
    
    pb_Qt_Z(i, :) = charge_poisson(pb_X, pb_R, P0_values);
    mpb_Qt_Z(i, :) = charge_poisson(mpb_X, mpb_R, P0_values);
end

pb_C_Z = zeros(size(pb_Qt_Z));
mpb_C_Z = zeros(size(mpb_Qt_Z));

for i = 1 : 2
    pb_C_Z(i, :) = pb_Qt_Z(i, :) ./ P0_values;
    mpb_C_Z(i, :) = mpb_Qt_Z(i, :) ./ P0_values;
end

dV = P0_values(2) - P0_values(1);
pb_Cdiff_Z = zeros(size(pb_Qt_Z));
mpb_Cdiff_Z = zeros(size(mpb_Qt_Z));

for i = 1 : 2
    for j = 1 : numel(P0_values) - 1
        pb_Cdiff_Z(i, j + 1) = (pb_Qt_Z(i, j + 1) - pb_Qt_Z(i, j)) / dV;
        mpb_Cdiff_Z(i, j + 1) = (mpb_Qt_Z(i, j + 1) - mpb_Qt_Z(i, j)) / dV;
    end
end

