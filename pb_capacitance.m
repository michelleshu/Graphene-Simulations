% PB_3P_CAPACITANCE Compute capacitance v. applied potential function from
% total charge function.

% Capacitance: C = Q / V

[P0, Q] = pb_3p_charge_poisson;
C = Q ./ P0;
plot(P0, C);
title(['Capacitance v. Applied Potential (P0)',...
    '(z = 1, C0 = 0.1 M, er = 80)'], 'FontSize', 16);
xlabel('P0 (V)', 'FontSize', 16);
ylabel('Capacitance (F)', 'FontSize', 16);
pause;

% Differential Capacitance: Cdiff = dQ / dV;

dV = P0(2) - P0(1);
Cdiff = zeros(size(Q));
for i = 1 : numel(Q) - 1
    Cdiff(i + 1) = (Q(i + 1) - Q(i)) / dV;
end
plot(P0(1, 2 : end), Cdiff(1, 2 : end));
title(['Differential Capacitance v. Applied Potential (P0)',...
    '(z = 1, C0 = 0.1 M, er = 80)'], 'FontSize', 16);
xlabel('P0 (V)', 'FontSize', 16);
ylabel('Differential Capacitance (F)', 'FontSize', 16);

pause;
semilogy(P0(1, 2 : end), Cdiff(1, 2 : end));