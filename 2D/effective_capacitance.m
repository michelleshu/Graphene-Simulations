% Compute effective capacitance of entire graphene surface. First,
% calculate capacitance value at each individual point along channel and
% then average the results.
function C_eff = effective_capacitance(X, P, R)

    E_0 = 8.854187817e-12;  % Vacuum permittivity (F/m)  
    E_R = 78.3;             % Relative permittivity (80 for H2O)

    C = zeros(size(R, 1), 1);

    for i = 2 : size(R, 1)
        % Approximate integral by trapezoidal estimation for total charge Qt
        C(i) = trapz(X, R(i, :) .* (E_R * E_0)) / P(i);
    end
    C_eff = mean(C);  % convert to F/m^2
    
end