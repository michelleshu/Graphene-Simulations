% MPB_CHARGE_CAP_PARAMS.M
% Compute functions for charge and capacitance due to Modified Poisson
% Boltzmann model while varying the following parameters: 
% bulk concentration (0.01 - 0.3 M), relative permittivity (30 - 120),
% effective ion size (0.1 - 5 nm), valence (+/- 1 or +/- 2)

% Michelle Shu | June 23, 2013

P0_values = (0.01 : 0.02 : 1.01);

% Defaults
C_0 = 150;  % mol/m^3
E_R = 78.3;
EFF = 1e-9;
Z = 1;

% Options for parameter variance
C_0_options = [10 50 100 150 300]; % in mol/m^3 = M * 1e3
E_R_options = [30 50 78.3 100 120];
EFF_options = [1e-10 1e-9 2e-9 3e-9 5e-9];
Z_options = [1 2];

% A. Bulk Concentration

Qt_C_0 = zeros(numel(C_0_options), numel(P0_values));

for i = 1 : numel(C_0_options)
    P = zeros(numel(P0_values), 301);
    R = zeros(numel(P0_values), 301);
    
    for j = 1 : numel(P0_values)
        disp([num2str(i) '   ' num2str(j)]);
        [X, P(j, :), R(j, :)] = mpb_potential(P0_values(j),...
            C_0_options(i), E_R, EFF, Z);
    end
    
    Qt_C_0(i, :) = charge_poisson(X, R, P0_values);
end

semilogy(P0_values, Qt_C_0(1, :) ./ 1e4, '-r', ...
         P0_values, Qt_C_0(2, :) ./ 1e4, '-g', ...
         P0_values, Qt_C_0(3, :) ./ 1e4, '-b', ...
         P0_values, Qt_C_0(4, :) ./ 1e4, '-c', ...
         P0_values, Qt_C_0(5, :) ./ 1e4, '-m');
title(['Total Charge per Area v. Applied Voltage (P0)',...
    ' (z = 1, er = 78.3, ion size = 1 nm)'], 'FontSize', 16);
xlabel('P0 (V)', 'FontSize', 16);
ylabel('Charge per Area (C / cm^2)', 'FontSize', 16);
legend('C0 = 0.01 M', 'C0 = 0.05 M', 'C0 = 0.10 M', 'C0 = 0.15 M',...
    'C0 = 0.30 M');

pause;

C_C_0 = zeros(size(Qt_C_0));

for i = 1 : 5
    C_C_0(i, :) = Qt_C_0(i, :) ./ P0_values;
end

semilogy(P0_values, C_C_0(1, :) ./ 1e4, '-r', ...
         P0_values, C_C_0(2, :) ./ 1e4, '-g', ...
         P0_values, C_C_0(3, :) ./ 1e4, '-b', ...
         P0_values, C_C_0(4, :) ./ 1e4, '-c', ...
         P0_values, C_C_0(5, :) ./ 1e4, '-m');
title(['Capacitance per Area v. Applied Potential (P0) ',...
    '(z = 1, er = 78.3, ion size = 1 nm)'], 'FontSize', 16);
xlabel('P0 (V)', 'FontSize', 16);
ylabel('Capacitance per Area (F / cm^2)', 'FontSize', 16);
legend('C0 = 0.01 M', 'C0 = 0.05 M', 'C0 = 0.10 M', 'C0 = 0.15 M',...
    'C0 = 0.30 M');

pause;

dV = P0_values(2) - P0_values(1);
Cdiff_C_0 = zeros(size(Qt_C_0));

for i = 1 : 5
    for j = 1 : numel(P0_values) - 1
        Cdiff_C_0(i, j + 1) = (Qt_C_0(i, j + 1) - Qt_C_0(i, j)) / dV;
    end
end

semilogy(P0_values(1, 2 : end), Cdiff_C_0(1, 2 : end) ./ 1e4, '-r', ...
         P0_values(1, 2 : end), Cdiff_C_0(2, 2 : end) ./ 1e4, '-g', ...
         P0_values(1, 2 : end), Cdiff_C_0(3, 2 : end) ./ 1e4, '-b', ...
         P0_values(1, 2 : end), Cdiff_C_0(4, 2 : end) ./ 1e4, '-c', ...
         P0_values(1, 2 : end), Cdiff_C_0(5, 2 : end) ./ 1e4, '-m');
        
title(['Differential Capacitance per Area v. Applied Potential (P0) ',...
    '(z = 1, er = 78.3, ion size = 1 nm)'], 'FontSize', 16);
xlabel('P0 (V)', 'FontSize', 16);
ylabel('Differential Capacitance (F / cm ^2)', 'FontSize', 16);
legend('C0 = 0.01 M', 'C0 = 0.05 M', 'C0 = 0.10 M', 'C0 = 0.15 M',...
    'C0 = 0.30 M');

% B. Relative Permittivity

Qt_E_R = zeros(numel(E_R_options), numel(P0_values));

for i = 1 : numel(E_R_options)
    P = zeros(numel(P0_values), 301);
    R = zeros(numel(P0_values), 301);
    
    for j = 1 : numel(P0_values)
        disp([num2str(i) '   ' num2str(j)]);
        [X, P(j, :), R(j, :)] = mpb_potential(P0_values(j),...
            C_0, E_R_options(i), EFF, Z);
    end
    
    Qt_E_R(i, :) = charge_poisson(X, R, P0_values);
end

semilogy(P0_values, Qt_E_R(1, :) ./ 1e4, '-r', ...
         P0_values, Qt_E_R(2, :) ./ 1e4, '-g', ...
         P0_values, Qt_E_R(3, :) ./ 1e4, '-b', ...
         P0_values, Qt_E_R(4, :) ./ 1e4, '-c', ...
         P0_values, Qt_E_R(5, :) ./ 1e4, '-m');
title(['Total Charge per Area v. Applied Voltage (P0)',...
    ' (z = 1, C0 = 0.15 M, ion size = 1 nm)'], 'FontSize', 16);
xlabel('P0 (V)', 'FontSize', 16);
ylabel('Charge per Area (C / cm^2)', 'FontSize', 16);
legend('er = 30', 'er = 50', 'er = 78.3', 'er = 100', 'er = 120');

pause;

C_E_R = zeros(size(Qt_E_R));

for i = 1 : 5
    C_E_R(i, :) = Qt_E_R(i, :) ./ P0_values;
end

semilogy(P0_values, C_E_R(1, :) ./ 1e4, '-r', ...
         P0_values, C_E_R(2, :) ./ 1e4, '-g', ...
         P0_values, C_E_R(3, :) ./ 1e4, '-b', ...
         P0_values, C_E_R(4, :) ./ 1e4, '-c', ...
         P0_values, C_E_R(5, :) ./ 1e4, '-m');
title(['Capacitance per Area v. Applied Potential (P0) ',...
    '(z = 1, C0 = 0.15 M, ion size = 1 nm)'], 'FontSize', 16);
xlabel('P0 (V)', 'FontSize', 16);
ylabel('Capacitance per Area (F / cm^2)', 'FontSize', 16);
legend('er = 30', 'er = 50', 'er = 78.3', 'er = 100', 'er = 120');

pause;

dV = P0_values(2) - P0_values(1);
Cdiff_E_R = zeros(size(Qt_E_R));

for i = 1 : 5
    for j = 1 : numel(P0_values) - 1
        Cdiff_E_R(i, j + 1) = (Qt_E_R(i, j + 1) - Qt_E_R(i, j)) / dV;
    end
end

semilogy(P0_values(1, 2 : end), Cdiff_E_R(1, 2 : end) ./ 1e4, '-r', ...
         P0_values(1, 2 : end), Cdiff_E_R(2, 2 : end) ./ 1e4, '-g', ...
         P0_values(1, 2 : end), Cdiff_E_R(3, 2 : end) ./ 1e4, '-b', ...
         P0_values(1, 2 : end), Cdiff_E_R(4, 2 : end) ./ 1e4, '-c', ...
         P0_values(1, 2 : end), Cdiff_E_R(5, 2 : end) ./ 1e4, '-m');
        
title(['Differential Capacitance per Area v. Applied Potential (P0) ',...
    '(z = 1, C0 = 0.15 M, ion size = 1 nm)'], 'FontSize', 16);
xlabel('P0 (V)', 'FontSize', 16);
ylabel('Differential Capacitance (F / cm ^2)', 'FontSize', 16);
legend('er = 30', 'er = 50', 'er = 78.3', 'er = 100', 'er = 120');

% C. Effective Ion Size

Qt_EFF = zeros(numel(EFF_options), numel(P0_values));

for i = 1 : numel(EFF_options)
    P = zeros(numel(P0_values), 301);
    R = zeros(numel(P0_values), 301);
    
    for j = 1 : numel(P0_values)
        disp([num2str(i) '   ' num2str(j)]);
        [X, P(j, :), R(j, :)] = mpb_potential(P0_values(j),...
            C_0, E_R, EFF_options(i), Z);
    end
    
    Qt_EFF(i, :) = charge_poisson(X, R, P0_values);
end

semilogy(P0_values, Qt_EFF(1, :) ./ 1e4, '-r', ...
         P0_values, Qt_EFF(2, :) ./ 1e4, '-g', ...
         P0_values, Qt_EFF(3, :) ./ 1e4, '-b', ...
         P0_values, Qt_EFF(4, :) ./ 1e4, '-c', ...
         P0_values, Qt_EFF(5, :) ./ 1e4, '-m');
title(['Total Charge per Area v. Applied Voltage (P0)',...
    ' (z = 1, C0 = 0.15 M, er = 78.3)'], 'FontSize', 16);
xlabel('P0 (V)', 'FontSize', 16);
ylabel('Charge per Area (C / cm^2)', 'FontSize', 16);
legend('ion size = 0.1 nm', 'ion size = 1.0 nm', 'ion size = 2.0 nm',...
    'ion size = 3.0 nm', 'ion size = 5.0 nm');
pause;

C_EFF = zeros(size(Qt_EFF));

for i = 1 : 5
    C_EFF(i, :) = Qt_EFF(i, :) ./ P0_values;
end

semilogy(P0_values, C_EFF(1, :) ./ 1e4, '-r', ...
         P0_values, C_EFF(2, :) ./ 1e4, '-g', ...
         P0_values, C_EFF(3, :) ./ 1e4, '-b', ...
         P0_values, C_EFF(4, :) ./ 1e4, '-c', ...
         P0_values, C_EFF(5, :) ./ 1e4, '-m');
title(['Capacitance per Area v. Applied Potential (P0) ',...
    '(z = 1, C0 = 0.15 M, er = 78.3)'], 'FontSize', 16);
xlabel('P0 (V)', 'FontSize', 16);
ylabel('Capacitance per Area (F / cm^2)', 'FontSize', 16);
legend('ion size = 0.1 nm', 'ion size = 1.0 nm', 'ion size = 2.0 nm',...
    'ion size = 3.0 nm', 'ion size = 5.0 nm');

pause;

dV = P0_values(2) - P0_values(1);
Cdiff_EFF = zeros(size(Qt_EFF));

for i = 1 : 5
    for j = 1 : numel(P0_values) - 1
        Cdiff_EFF(i, j + 1) = (Qt_EFF(i, j + 1) - Qt_EFF(i, j)) / dV;
    end
end

semilogy(P0_values(1, 2 : end), Cdiff_EFF(1, 2 : end) ./ 1e4, '-r', ...
         P0_values(1, 2 : end), Cdiff_EFF(2, 2 : end) ./ 1e4, '-g', ...
         P0_values(1, 2 : end), Cdiff_EFF(3, 2 : end) ./ 1e4, '-b', ...
         P0_values(1, 2 : end), Cdiff_EFF(4, 2 : end) ./ 1e4, '-c', ...
         P0_values(1, 2 : end), Cdiff_EFF(5, 2 : end) ./ 1e4, '-m');
        
title(['Differential Capacitance per Area v. Applied Potential (P0) ',...
    '(z = 1, C0 = 0.15 M, er = 78.3)'], 'FontSize', 16);
xlabel('P0 (V)', 'FontSize', 16);
ylabel('Differential Capacitance (F / cm ^2)', 'FontSize', 16);
legend('ion size = 0.1 nm', 'ion size = 1.0 nm', 'ion size = 2.0 nm',...
    'ion size = 3.0 nm', 'ion size = 5.0 nm');

% D. Valence

Qt_Z = zeros(numel(Z_options), numel(P0_values));
 
for i = 1 : numel(Z_options)
    P = zeros(numel(P0_values), 301);
    R = zeros(numel(P0_values), 301);
    
    for j = 1 : numel(P0_values)
        disp([num2str(i) '   ' num2str(j)]);
        [X, P(j, :), R(j, :)] = mpb_potential(P0_values(j),...
            C_0, E_R, EFF, Z_options(i));
    end
    
    Qt_Z(i, :) = charge_poisson(X, R, P0_values);
end

semilogy(P0_values, Qt_Z(1, :) ./ 1e4, '-r', ...
         P0_values, Qt_Z(2, :) ./ 1e4, '-b');
title(['Total Charge per Area v. Applied Voltage (P0) ',...
    '(C0 = 0.15 M, er = 78.3, ion size = 1 nm)'], 'FontSize', 16);
xlabel('P0 (V)', 'FontSize', 16);
ylabel('Charge per Area (C / cm^2)', 'FontSize', 16);
legend('Z = 1', 'Z = 2');
pause;

C_Z = zeros(size(Qt_Z));

for i = 1 : 2
    C_Z(i, :) = Qt_Z(i, :) ./ P0_values;
end

semilogy(P0_values, C_Z(1, :) ./ 1e4, '-r', ...
         P0_values, C_Z(2, :) ./ 1e4, '-b');
title(['Capacitance per Area v. Applied Potential (P0) ',...
    '(C0 = 0.15 M, er = 78.3, ion size = 1 nm)'], 'FontSize', 16);
xlabel('P0 (V)', 'FontSize', 16);
ylabel('Capacitance per Area (F / cm^2)', 'FontSize', 16);
legend('Z = 1', 'Z = 2');

pause;

dV = P0_values(2) - P0_values(1);
Cdiff_Z = zeros(size(Qt_Z));

for i = 1 : 2
    for j = 1 : numel(P0_values) - 1
        Cdiff_Z(i, j + 1) = (Qt_Z(i, j + 1) - Qt_Z(i, j)) / dV;
    end
end

semilogy(P0_values(1, 2 : end), Cdiff_Z(1, 2 : end) ./ 1e4, '-r', ...
         P0_values(1, 2 : end), Cdiff_Z(2, 2 : end) ./ 1e4, '-b');
        
title(['Differential Capacitance per Area v. Applied Potential (P0) ',...
    '(C0 = 0.15 M, er = 78.3, ion size = 1 nm)'], 'FontSize', 16);
xlabel('P0 (V)', 'FontSize', 16);
ylabel('Differential Capacitance (F / cm ^2)', 'FontSize', 16);
legend('Z = 1', 'Z = 2');

