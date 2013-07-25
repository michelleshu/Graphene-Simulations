function C = C_top( V )
% C_TOP Exponential fit from MATLAB of capacitance as a function of 
% applied voltage.  
    C = 0.546 * exp(-16.65 .* V) + 0.2203 * exp(-0.9717 .* V);
end

