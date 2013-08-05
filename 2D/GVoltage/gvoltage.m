% GVOLTAGE Find function of potential along graphene surface using
% MATLAB's fsolve function
% Michelle Shu | July 24, 2013

function [V, fval, X] = gvoltage(W, L, VGS_TOP, VSD)

    I = current(W, L, VGS_TOP, VSD);
    H = 1e-7;           % Resolution (grid size along channel length)

    % Compute potential function
    X = (0 : H : L);
    V_init = (VSD / L) .* X; % Initial guess for V

    f = @(V)V_functionsB(V, I, W, H, VGS_TOP, VSD);

    options = optimoptions('fsolve', 'Algorithm', 'levenberg-marquardt', ...
        'Display', 'iter', 'TolFun', 1e-5, 'TolX', 1e-5);
    [V, fval] = fsolve(f, V_init, options);

end