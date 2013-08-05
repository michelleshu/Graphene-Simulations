function C = C_top_3nmEFF( V )
% C_TOP_3nmEFF Exponential fit from MATLAB of capacitance as a function of 
% applied voltage. Effective ion size = 3 nm

% This is an exponential curve fit, accurate between 0 and 1 V.
    C = 0.546 * exp(-16.65 .* V) + 0.2203 * exp(-0.9717 .* V);
end