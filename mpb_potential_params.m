% MPB_POTENTIAL_PARAMS Script
% Modified Poisson Boltzmann Model
% Display how electric potential profile varies with: 
% bulk concentration (0.01 - 0.3 M), relative permittivity (30 - 120),
% effective ion size (0.1 - 5 nm)

% Michelle Shu | June 23, 2013

% A. Bulk Concentration

C_0_options = [10 50 100 150 300]; % in mol/m^3 = M * 1e3
P_0 = 0.025;
E_R = 78.3;
EFF = 3e-10;
% 
% [X1, P1, ~] = mpb_potential(P_0, C_0_options(1), E_R, EFF);
% [X2, P2, ~] = mpb_potential(P_0, C_0_options(2), E_R, EFF);
% [X3, P3, ~] = mpb_potential(P_0, C_0_options(3), E_R, EFF);
% [X4, P4, ~] = mpb_potential(P_0, C_0_options(4), E_R, EFF);
% [X5, P5, ~] = mpb_potential(P_0, C_0_options(5), E_R, EFF);
% 
% plot(X1, P1, '-r', X2, P2, '-g', X3, P3, '-b', X4, P4, '-c', X5, P5, '-m');
% title(['Potential v. Distance from Interface (MPB) ',...
%     '(z = 1, P0 = 25 mV, er = 80, ion size = 0.3 nm)'], 'FontSize', 16);
% xlabel('Distance (m)', 'FontSize', 16);
% ylabel('Potential (V)', 'FontSize', 16);
% legend('C0 = 0.01 M', 'C0 = 0.05 M', 'C0 = 0.10 M', 'C0 = 0.15 M',...
%     'C0 = 0.30 M');
% 
% pause;

% B. Permittivity

C_0 = 150;
E_R_options = [30 50 78.3 100 120];

[X1, P1, ~] = mpb_potential(P_0, C_0, E_R_options(1), EFF);
[X2, P2, ~] = mpb_potential(P_0, C_0, E_R_options(2), EFF);
[X3, P3, ~] = mpb_potential(P_0, C_0, E_R_options(3), EFF);
[X4, P4, ~] = mpb_potential(P_0, C_0, E_R_options(4), EFF);
[X5, P5, ~] = mpb_potential(P_0, C_0, E_R_options(5), EFF);


plot(X1, P1, '-r', X2, P2, '-g', X3, P3, '-b', X4, P4, '-c', X5, P5, '-m');
title(['Potential v. Distance from Interface (MPB) ',...
    '(z = 1, P0 = 25 mV, C0 = 0.15 M, ion size = 0.3 nm)'], 'FontSize', 16);
xlabel('Distance (m)', 'FontSize', 16);
ylabel('Potential (V)', 'FontSize', 16);
legend('er = 30', 'er = 50', 'er = 78.3', 'er = 100', 'er = 120');

pause;

% C. Effective Ion Size

EFF_options = [1e-10 1e-9 2e-9 3e-9 5e-9];

[X1, P1, ~] = mpb_potential(P_0, C_0, E_R, EFF_options(1));
[X2, P2, ~] = mpb_potential(P_0, C_0, E_R, EFF_options(2));
[X3, P3, ~] = mpb_potential(P_0, C_0, E_R, EFF_options(3));
[X4, P4, ~] = mpb_potential(P_0, C_0, E_R, EFF_options(4));
[X5, P5, ~] = mpb_potential(P_0, C_0, E_R, EFF_options(5));

plot(X1, P1, '-r', X2, P2, '-g', X3, P3, '-b', X4, P4, '-c', X5, P5, '-m');
title(['Potential v. Distance from Interface (MPB) ',...
    '(z = 1, C0 = 0.15 M, P0 = 25 mV, er = 80)'], 'FontSize', 16);
xlabel('Distance (m)', 'FontSize', 16);
ylabel('Potential (V)', 'FontSize', 16);
legend('ion size = 0.1 nm', 'ion size = 1.0 nm', 'ion size = 2.0 nm', ...
    'ion size = 3.0 nm', 'ion size = 5.0 nm');


