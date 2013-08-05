% GVOLTAGE Find function of potential along graphene surface using
% MATLAB's fsolve function
% Michelle Shu | July 24, 2013

function [V, fval, X] = gvoltage(W, L, V_0, VGS_TOP)

    I = current(W, L, VGS_TOP);
    H = 1e-7; % Resolution (grid size along channel length)

    % Compute potential function
    X = (0 : H : L);
    V_init = (V_0 / L) .* X; % Initial guess for V

    f = @(V)V_functions(V, I, W, H, V_0, VGS_TOP);

    options = optimoptions('fsolve', 'Algorithm', 'levenberg-marquardt', ...
        'Display', 'iter', 'TolFun', 1e-6, 'TolX', 1e-6);
    [V, fval] = fsolve(f, V_init, options);

end