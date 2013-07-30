% GVOLTAGE Find function of potential along graphene surface using
% MATLAB's fsolve function
% Michelle Shu | July 24, 2013

H = 5e-9;
L = 3e-5;

% PART 2: Compute potential function
V_0 = 0.1;               % Applied potential
X = (0 : H : L);
V_init = (V_0 / L) .* X; % Initial guess for V

options = optimoptions('fsolve', 'Algorithm', 'levenberg-marquardt', ...
    'Display', 'iter', 'TolFun', 1e-6, 'TolX', 1e-6);
[V, fval] = fsolve(@V_functions, V_init, options);    