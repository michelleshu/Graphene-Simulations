% GVOLTAGE Find function of potential along graphene surface
% Michelle Shu | July 22, 2013

% PART 1: Compute current I
W = 5e-6;           % Channel width
L = 3e-5;           % Sensor length
H = 1e-9;           % Grid size along sensor length

Q = 1.60217657e-19; % Elementary charge
N_0 = 5e15;         % Minimum carrier density
MOB = 0.055;        % m^2 / s

VGS_TOP = 0.4;      % Volts
K = 1e-3;           % Grid size along V
V_values = (0 : K : VGS_TOP); 

% C_top = polyval(polyFitCapacitance, V_values);   % polynomial fit
C_top = feval(expFitCapacitance, V_values)';    % exponential fit

F1 = (C_top .* (VGS_TOP - V_values) / Q) .^ 2;
F2 = sqrt((N_0^2) + F1);
I = - Q * W * MOB/ L * trapz(F2, V_values); % Derive I by integration

% PART 2: Compute potential function
V_0 = 0.1;            % Applied potential
X = (0 : H : L);
JUMP = 1e-5;          % Jump for relaxation
CONV = 0.01;          % Max difference for convergencex

% Initialization
V = (V_0 / L) .* X;
done = false;

% Loop until convergence
while ~done
    max_diff = 0;
    for i = 2 : numel(V) - 1
        % Compute dV/dx from central derivative theorem
        dVdx = (V(i + 1) - V(i - 1)) / (2 * H);
        % Generate function to solve for V from
        F3 = I - (Q * W * MOB * dVdx * F2);   % equals 0
        % Use bisection method to solve for V
        rindex = findrootindex(F3);
        if (rindex == -1) % returns -1 if no root found
              plot(V_values, F4);
              pause;
              plot(X, V);
              pause;
        end;
        V_calc = V_values(rindex);
        disp([V_calc V(i)]);
        % Record diff if max
        if (abs(V_calc - V(i)) / V(i)) > max_diff
            max_diff = abs(V_calc - V(i) / V(i));
        end
        
        % Relax
        V(i) = V(i) + JUMP * (V_calc - V(i));
    end
    % disp(max_diff);
    % Convergence check
    if (max_diff <= CONV)
        done = true;
    end
end