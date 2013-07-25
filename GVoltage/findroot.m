% Find root via bisection method.
% Michelle Shu | July 19, 2013

function root = findroot(F)
    min = 0;
    max = numel(F);
    while true
        % Base cases
        if (min == max)
            root = min;
            return;   
        elseif (max - min == 1)
            root = (max + min) / 2;
            return;
        else
            mid = (min + max) / 2;
            if ((F(min) * F(mid)) < 0)
                max = mid;
            elseif ((F(mid) * F(max)) < 0)
                min = mid;
            else
                root = -1;  % Error: No zeros
                return;
            end
        end
    end
end