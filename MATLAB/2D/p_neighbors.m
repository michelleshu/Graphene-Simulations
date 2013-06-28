function p_sum = p_neighbors ( P, x, y, x_max, y_max )
% GET_P_NEIGHBORS Retrieve potential values in matrix P of the 4 neighbors
% of entry (x, y). Return the sum

% Michelle Shu | June 24, 2013

% Left Neighbor
if (x > 1)
    p_left = P(y, x - 1);
else
    p_left = P(y, x_max);
end

% Right Neighbor
if (x < x_max)
    p_right = P(y, x + 1);
else
    p_right = P(y, 1);
end

% Top Neighbor
if (y > 1)
    p_top = P(y - 1, x);
else
    p_top = P(y_max, x);
end

% Bottom Neighbor
if (y < y_max)
    p_bottom = P(y + 1, x);
else
    p_bottom = P(1, x);
end

p_sum = p_left + p_right + p_top + p_bottom;

end

