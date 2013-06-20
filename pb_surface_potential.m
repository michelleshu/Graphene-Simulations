P_0_options = [0.025, 0.05, 0.075, 0.1];
C_0 = 100;
E_R = 80;

[X1, P1] = pb_3p_double_layer(P_0_options(1), C_0, E_R);
[X2, P2] = pb_3p_double_layer(P_0_options(2), C_0, E_R);
[X3, P3] = pb_3p_double_layer(P_0_options(3), C_0, E_R);
[X4, P4] = pb_3p_double_layer(P_0_options(4), C_0, E_R);

plot(X1, P1, '-r', X2, P2, '-g', X3, P3, '-b', X4, P4, '-c');
title(['Potential v. Distance from Interface',...
    '(z = 1, C0 = 0.1 M, er = 80)'], 'FontSize', 16);
xlabel('Distance (m)', 'FontSize', 16);
ylabel('Potential (V)', 'FontSize', 16);
legend('P0 = 25 mV', 'P0 = 50 mV', 'P0 = 75 mV', 'P0 = 100mV');