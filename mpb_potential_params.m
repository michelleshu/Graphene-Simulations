% MPB_POTENTIAL_PARAMS Script
% Modified Poisson Boltzmann Model
% Display how electric potential profile varies with: 
% bulk concentration (0.01 - 0.3 M), relative permittivity (30 - 120),
% effective ion size (0.1 - 5 nm), valence (+/- 1 or +/- 2)

% Michelle Shu | June 23, 2013

% Defaults
P_0 = 0.025;
C_0 = 150;
E_R = 78.3;
EFF = 1e-9;
Z = 1;

% Options for parameter variance
C_0_options = [10 50 100 150 300]; % in mol/m^3 = M * 1e3
E_R_options = [30 50 78.3 100 120];
EFF_options = [1e-10 1e-9 2e-9 3e-9 5e-9];
Z_options = [1 2];

% A. Bulk Concentration

[X1, P1, ~] = mpb_potential(P_0, C_0_options(1), E_R, EFF, Z);
[X2, P2, ~] = mpb_potential(P_0, C_0_options(2), E_R, EFF, Z);
[X3, P3, ~] = mpb_potential(P_0, C_0_options(3), E_R, EFF, Z);
[X4, P4, ~] = mpb_potential(P_0, C_0_options(4), E_R, EFF, Z);
[X5, P5, ~] = mpb_potential(P_0, C_0_options(5), E_R, EFF, Z);

plot(X1, P1, '-r', X2, P2, '-g', X3, P3, '-b', X4, P4, '-c', X5, P5, '-m');
title(['Potential v. Distance from Interface (MPB) ',...
    '(z = 1, P0 = 25 mV, er = 80, ion size = 1 nm)'], 'FontSize', 16);
xlabel('Distance (m)', 'FontSize', 16);
ylabel('Potential (V)', 'FontSize', 16);
legend('C0 = 0.01 M', 'C0 = 0.05 M', 'C0 = 0.10 M', 'C0 = 0.15 M',...
    'C0 = 0.30 M');

pause;

% B. Permittivity

[X1, P1, ~] = mpb_potential(P_0, C_0, E_R_options(1), EFF, Z);
[X2, P2, ~] = mpb_potential(P_0, C_0, E_R_options(2), EFF, Z);
[X3, P3, ~] = mpb_potential(P_0, C_0, E_R_options(3), EFF, Z);
[X4, P4, ~] = mpb_potential(P_0, C_0, E_R_options(4), EFF, Z);
[X5, P5, ~] = mpb_potential(P_0, C_0, E_R_options(5), EFF, Z);


plot(X1, P1, '-r', X2, P2, '-g', X3, P3, '-b', X4, P4, '-c', X5, P5, '-m');
title(['Potential v. Distance from Interface (MPB) ',...
    '(z = 1, P0 = 25 mV, C0 = 0.15 M, ion size = 1 nm)'], 'FontSize', 16);
xlabel('Distance (m)', 'FontSize', 16);
ylabel('Potential (V)', 'FontSize', 16);
legend('er = 30', 'er = 50', 'er = 78.3', 'er = 100', 'er = 120');

pause;

% C. Effective Ion Size

[X1, P1, ~] = mpb_potential(P_0, C_0, E_R, EFF_options(1), Z);
[X2, P2, ~] = mpb_potential(P_0, C_0, E_R, EFF_options(2), Z);
[X3, P3, ~] = mpb_potential(P_0, C_0, E_R, EFF_options(3), Z);
[X4, P4, ~] = mpb_potential(P_0, C_0, E_R, EFF_options(4), Z);
[X5, P5, ~] = mpb_potential(P_0, C_0, E_R, EFF_options(5), Z);

plot(X1, P1, '-r', X2, P2, '-g', X3, P3, '-b', X4, P4, '-c', X5, P5, '-m');
title(['Potential v. Distance from Interface (MPB) ',...
    '(z = 1, C0 = 0.15 M, P0 = 25 mV, er = 80)'], 'FontSize', 16);
xlabel('Distance (m)', 'FontSize', 16);
ylabel('Potential (V)', 'FontSize', 16);
legend('ion size = 0.1 nm', 'ion size = 1.0 nm', 'ion size = 2.0 nm', ...
    'ion size = 3.0 nm', 'ion size = 5.0 nm');

pause;

% D. Valence

[X1, P1, ~] = mpb_potential(P_0, C_0, E_R, EFF, Z_options(1));
[X2, P2, ~] = mpb_potential(P_0, C_0, E_R, EFF, Z_options(2));

plot(X1, P1, '-r', X2, P2, '-b');
title(['Potential v. Distance from Interface (MPB) ',...
    '(C0 = 0.15 M, P0 = 25 mV, er = 80, ion size = 1 nm)'], 'FontSize', 16);
xlabel('Distance (m)', 'FontSize', 16);
ylabel('Potential (V)', 'FontSize', 16);
legend('Z = 1', 'Z = 2');


