function C_eff = CV_curve()

    mex potential.c;

    P0_values = [0.1 : 0.1 : 1];    % V
    C_eff = zeros(size(P0_values));
    E_R = 78.3;
    EFF = 3e-9;
    CONC = 0.1;

    for i = 1 : numel(P0_values)
        [X, Y, P, R] = potential(P0_values(i), E_R, EFF, 1, 1, CONC, 1);
        C_eff(i) = effective_capacitance(X, P, R);
        disp(i);
    end
    
end