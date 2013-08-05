function C = C_top_5nmEFF( V )
% C_TOP_3nmEFF Exponential fit from MATLAB of capacitance as a function of 
% applied voltage. Effective ion size = 3 nm

% This is an exponential curve fit, accurate between 0 and 1 V.
    C = 0.7389 * exp(-40.22 * V) + 0.1869 * exp(-2.394 * V);
end