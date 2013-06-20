C_0_options = [0.1, 1, 10, 100, 1000]; % in mol/m^3
P_0 = 0.025;
E_R = 80;

[X1, P1] = pb_3p_double_layer(P_0, C_0_options(1), E_R);
[X2, P2] = pb_3p_double_layer(P_0, C_0_options(2), E_R);
[X3, P3] = pb_3p_double_layer(P_0, C_0_options(3), E_R);
[X4, P4] = pb_3p_double_layer(P_0, C_0_options(4), E_R);
[X5, P5] = pb_3p_double_layer(P_0, C_0_options(5), E_R);

plot(X1, P1, '-r', X2, P2, '-g', X3, P3, '-b', X4, P4, '-c', X5, P5, '-m');
title(['Potential v. Distance from Interface',...
    '(z = 1, P0 = 25 mV, er = 80)'], 'FontSize', 16);
xlabel('Distance (m)', 'FontSize', 16);
ylabel('Potential (V)', 'FontSize', 16);
legend('C0 = 1e-4 M', 'C0 = 1e-3 M', 'C0 = 1e-2 M', 'C0 = 0.1 M', 'C0 = 1 M');