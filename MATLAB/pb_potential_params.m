% PB_POTENTIAL_PARAMS Script
% Display how electric potential profile varies with: bulk concentration,
% surface potential, permittivity

% Michelle Shu | June 18, 2013

% A. Bulk Concentration

C_0_options = [0.1, 1, 10, 100, 1000]; % in mol/m^3
P_0 = 0.025;
E_R = 78.3;

[X1, P1, ~] = pb_potential(P_0, C_0_options(1), E_R);
[X2, P2, ~] = pb_potential(P_0, C_0_options(2), E_R);
[X3, P3, ~] = pb_potential(P_0, C_0_options(3), E_R);
[X4, P4, ~] = pb_potential(P_0, C_0_options(4), E_R);
[X5, P5, ~] = pb_potential(P_0, C_0_options(5), E_R);

plot(X1, P1, '-r', X2, P2, '-g', X3, P3, '-b', X4, P4, '-c', X5, P5, '-m');
title(['Potential v. Distance from Interface',...
    '(z = 1, P0 = 25 mV, er = 80)'], 'FontSize', 16);
xlabel('Distance (m)', 'FontSize', 16);
ylabel('Potential (V)', 'FontSize', 16);
legend('C0 = 1e-4 M', 'C0 = 1e-3 M', 'C0 = 1e-2 M', 'C0 = 0.1 M', 'C0 = 1 M');

pause;

% B. Surface Potential

P_0_options = [0.025, 0.05, 0.075, 0.1];
C_0 = 100;
E_R = 80;

[X1, P1] = pb_potential(P_0_options(1), C_0, E_R);
[X2, P2] = pb_potential(P_0_options(2), C_0, E_R);
[X3, P3] = pb_potential(P_0_options(3), C_0, E_R);
[X4, P4] = pb_potential(P_0_options(4), C_0, E_R);

plot(X1, P1, '-r', X2, P2, '-g', X3, P3, '-b', X4, P4, '-c');
title(['Potential v. Distance from Interface',...
    '(z = 1, C0 = 0.1 M, er = 80)'], 'FontSize', 16);
xlabel('Distance (m)', 'FontSize', 16);
ylabel('Potential (V)', 'FontSize', 16);
legend('P0 = 25 mV', 'P0 = 50 mV', 'P0 = 75 mV', 'P0 = 100mV');

pause;

% C. Permittivity

E_R_options = [20 40 60 80 100];
P_0 = 0.025;
C_0 = 100;

[X1, P1] = pb_potential(P_0, C_0, E_R_options(1));
[X2, P2] = pb_potential(P_0, C_0, E_R_options(2));
[X3, P3] = pb_potential(P_0, C_0, E_R_options(3));
[X4, P4] = pb_potential(P_0, C_0, E_R_options(4));
[X5, P5] = pb_potential(P_0, C_0, E_R_options(5));

plot(X1, P1, '-r', X2, P2, '-g', X3, P3, '-b', X4, P4, '-c', X5, P5, '-m');
title(['Potential v. Distance from Interface',...
    '(z = 1, P0 = 25 mV, C0 = 0.1M)'], 'FontSize', 16);
xlabel('Distance (m)', 'FontSize', 16);
ylabel('Potential (V)', 'FontSize', 16);
legend('er = 20', 'er = 40', 'er = 60', 'er = 80', 'er = 100');