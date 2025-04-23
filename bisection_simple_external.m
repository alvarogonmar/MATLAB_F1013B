function output = bisection_simple_external(f, a, b, tol)


    % Check if the signs at the interval ends are opposite
    if f(a) * f(b) >= 0
        error('Function must have opposite signs at a and b');
    end

    counter = 0; % Initialize iteration counter
    
    while true
        m = (a + b) / 2; % Midpoint
        fm = f(m);
        
        % Check if the midpoint is a root or the error is within tolerance
        if fm == 0
            break;
        elseif f(a) * fm > 0
            a = m;
        else
            b = m;
        end

        counter = counter + 1;
        error_val = abs(b - a) / 2;
        if error_val < tol
            break;
        end
    end

    root = (a + b) / 2;
    output = [counter, root];
end
