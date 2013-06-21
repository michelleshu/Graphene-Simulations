% Constants
E_0 = 8.854187817e-12;  % Vacuum permittivity (F/m)
E_R = 78.3;

pb_Qt = zeros(1, 100);
mpb_Qt = zeros(1, 100);

for i = 1 : 100
    pb_q = pb_R(i, :) .* (E_R * E_0);
    mpb_q = mpb_R(i, :) .* (E_R * E_0);
    % Approximate integral by trapezoidal estimation for total charge Qt
    pb_Qt(i) = trapz(X, pb_q);
    mpb_Qt(i) = trapz(X, mpb_q);
end

semilogy(P0_values, pb_Qt, '-b', P0_values, mpb_Qt, '-g');
title(['Total Charge (Q) v. Applied Voltage (P0)',...
    '(z = 1, C0 = 0.1 M, er = 80)'], 'FontSize', 16);
xlabel('P0 (V)', 'FontSize', 16);
ylabel('Q (C)', 'FontSize', 16);
legend('PB', 'MPB');