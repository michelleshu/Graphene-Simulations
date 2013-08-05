function C = C_top_1nmEFF( V )
% C_TOP_1nmEFF Exponential fit from MATLAB of capacitance as a function of 
% applied voltage. Effective ion size = 1 nm

% This is a Gaussian curve fit, accurate between 0 and 1 V.
Ca = 0.6913 .* exp(-((V - 0.002266) ./ 0.3766) .^ 2);
Cb = -0.2307 .* exp(-((V + 0.007666) ./ 0.07477) .^ 2);
Cc = -0.2716 .* exp(-((V - 0.1425) ./ 0.2891) .^ 2);
Cd = 0.6821 .* exp(-((V + 0.5992) ./ 2.515) .^ 2);

C = Ca + Cb + Cc + Cd;

end