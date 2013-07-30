% NOTES: Fill in later to specify full routine

% Compute Vdrift
% Curve fit V0 and compute the derivative by hand to get E (electric field)
E = - (1135.8 * exp(55270 * X_5nm) + 813.28 * exp(-39270 * X_5nm));
Vdrift = 0.055 .* E;    % mobility * E
Vdrift_cm = Vdrift * 100;
plot(X_5nm .* 100, Vdrift_cm);

% Compute n0

