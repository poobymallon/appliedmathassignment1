function firstmethods()

    % initial guesses for each root-finding method
    bisect_start = [0, 5, 2.5];       % [left, right, midpoint] for bisection
    newton_start = 0;                 % Initial guess for Newton-Raphson
    secant_start = [0;1];             % Two initial points for secant method
    
    % Definition of the nonlinear test function
    test_func01 = @(x) (x.^3)/100 - (x.^2)/8 + 2*x + 6*sin(x/2+6) -.7 - exp(x/6);
    
    % Derivative of the test function for use in Newton-Raphson method
    test_derivative01 = @(x) 3*(x.^2)/100 - 2*x/8 + 2 +(6/2)*cos(x/2+6) - exp(x/6)/6;
    
    % Newton-Raphson update function
    newton_func = @(x) x - test_func01(x)/test_derivative01(x);
    
    % Secant method update function using two previous points
    secant_func = @(xi, xj) (xj*test_func01(xi)-xi*test_func01(xj))/(test_func01(xi)-test_func01(xj));
    
    % Apply each root-finding method
    [bi_root, bisect_list] = bisect_root_find(bisect_start, test_func01);
    [new_root, newton_list] = newton_root_find(newton_start, newton_func, test_func01);
    [sec_root, secant_list] = secant_root_find(secant_start, secant_func, test_func01);

    % Evaluate function over a range for plotting
    t = linspace(-15, 40, 2000);
    ft = test_func01(t);
    
    % Plot bisection method results
    figure
    subplot(1, 3, 1)
    hold on
    plot(t, ft)                              % Plot the function
    ms = bisect_list(:, 3);                  % Extract midpoints from iterations
    fms = test_func01(ms);                   % Function values at midpoints
    plot(ms, fms)                            % Plot midpoints
    plot(bi_root, test_func01(bi_root), "g*")% Final root estimate
    hold off
    
    % Plot Newton-Raphson method results
    subplot(1, 3, 2)
    hold on
    plot(t, ft)
    x1s = newton_list;                       % Iterative x-values
    f1s = test_func01(x1s);                  % Corresponding function values
    plot(x1s, f1s)
    plot(new_root, test_func01(new_root), "g*")
    hold off
    
    % Plot secant method results
    subplot(1, 3, 3)
    hold on
    plot(t, ft)
    x2s = secant_list;
    f2s = test_func01(x2s);
    plot(x2s, f2s)
    plot(sec_root, test_func01(sec_root), "g*")
    hold off

end

% One iteration of the bisection method to narrow the interval
function [Lnew, Rnew, mx] = bisect_func(input, func)
    L = input(1);
    R = input(2);
    if (func(L) > 0 && func(R) > 0) || (func(L) < 0 && func(R) < 0)
        return
    end
    mx = (L+R)/2;                  % Midpoint
    mf = func(mx);                 % Function value at midpoint
    Lf = func(L);                  % Function value at left
    Rf = func(R);                  % Function value at right

    % Determine the subinterval where the sign change occurs
    if Lf > 0 && mf < 0
        Lnew = L;
        Rnew = mx;
    end
    if Lf < 0 && mf > 0
        Lnew = L;
        Rnew = mx;
    end
    if Rf > 0 && mf < 0
        Lnew = mx;
        Rnew = R;
    end
    if Rf < 0 && mf > 0
        Lnew = mx;
        Rnew = R;
    end

end

% Complete bisection method until convergence
function [root, bisect_list] = bisect_root_find(input, func)
    bisect_list = zeros(1,3);               % Preallocate matrix for [L, R, midpoint]
    bisect_list(1,:) = input;
    mx = bisect_list(1, 3);
    iter = 1;

    % Iterate until the function value at midpoint is sufficiently small
    while abs(func(mx)) > 10^-14
        [L, R, mx] = bisect_func(bisect_list(iter, :), func);
        iter = iter+1;
        bisect_list(iter,:) = [L, R, mx];
    end
    root = mx;                              % Return the root estimate
end

% Newton-Raphson method until convergence
function [root, newton_list] = newton_root_find(start, newton_func, func)
    newton_list = zeros(1,1);               % Preallocate array for iterations
    newton_list(1) = start;
    newton_list(2) = newton_func(newton_list(1));
    iter = 2;

    % Iterate until convergence based on function value
    while abs(func(newton_list(iter))) > 10^-14
        iter = iter + 1;
        newton_list(iter) = newton_func(newton_list(iter-1));
    end
    root = newton_list(end);                % Final root estimate
end

% Secant method until convergence
function [root, secant_list] = secant_root_find(start, secant_func, func)
    secant_list = zeros(1,2);               % Preallocate for two starting values
    secant_list(1:2) = start;
    secant_list(3) = secant_func(secant_list(2), secant_list(1));
    iter = 3;

    % Iterate until convergence based on function value
    while abs(func(secant_list(iter))) > 10^-14
        iter = iter + 1;
        secant_list(iter) = secant_func(secant_list(iter-1), secant_list(iter-2));
    end
    root = secant_list(end);                % Final root estimate
end
