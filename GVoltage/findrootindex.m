% Find root via bisection method.
% Michelle Shu | July 19, 2013

function root = findrootindex(F)
    min = 1;
    max = numel(F);
    while true
        % Base cases
        if (min == max)
            root = min;
            return   
        elseif (max - min == 1)
            % Return the one with smallest absolute value.
            if abs(F(min)) < abs(F(max))
                root = min;
            else
                root = max;
            end
            return
        % Otherwise: Keep looking
        else
            mid = floor((min + max) / 2);
            if ((F(min) * F(mid)) <= 0)
                max = mid;
            elseif ((F(mid) * F(max)) <= 0)
                min = mid;
            else
                root = -1;  % Error: No zeros
                return
            end
        end
    end
end