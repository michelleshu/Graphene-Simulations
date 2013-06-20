E_R_options = [20 40 60 80 100];
P_0 = 0.025;
C_0 = 100;

[X1, P1] = pb_3p_double_layer(P_0, C_0, E_R_options(1));
[X2, P2] = pb_3p_double_layer(P_0, C_0, E_R_options(2));
[X3, P3] = pb_3p_double_layer(P_0, C_0, E_R_options(3));
[X4, P4] = pb_3p_double_layer(P_0, C_0, E_R_options(4));
[X5, P5] = pb_3p_double_layer(P_0, C_0, E_R_options(5));

plot(X1, P1, '-r', X2, P2, '-g', X3, P3, '-b', X4, P4, '-c', X5, P5, '-m');
title(['Potential v. Distance from Interface',...
    '(z = 1, P0 = 25 mV, C0 = 0.1M)'], 'FontSize', 16);
xlabel('Distance (m)', 'FontSize', 16);
ylabel('Potential (V)', 'FontSize', 16);
legend('er = 20', 'er = 40', 'er = 60', 'er = 80', 'er = 100');