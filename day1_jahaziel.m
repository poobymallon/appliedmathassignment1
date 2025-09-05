%% part 1
f = @(x) (x.^3)/100 - (x.^2)/8 + 2*x + 6*sin(x/2 + 6) - 0.7 - exp(x/6);
df = @(x) 3*(x.^2)/100 - x/4 + 2 + 3*cos(x/2 + 6) - exp(x/6)/6;

%where f(l) and f(r) have opposite signs
L = 0; R = 1;

% using this to check the signs are diff
fprintf('f(L) = %.4f, f(R) = %.4f\n', f(L), f(R));

figure;
fplot(f, [-15 40]);     
yline(0,'--');           % horizontal line at y=0
grid on;
title('checking where roots are');
xlabel('x'); ylabel('f(x)');

% midpoint for newton start
M = (L + R)/2;

% running all mthds
root_bisect = bisection_method(f, L, R, 1e-6);
root_newton = newton_method(f, df, M, 1e-6);
root_secant = secant_method(f, L, R, 1e-6);

% printing results
fprintf('bisection root: %.6f\n', root_bisect);
fprintf('newton root: %.6f\n', root_newton);
fprintf('secant root: %.6f\n', root_secant);

function root = bisection_method(f, L, R, tol)
    % usng interval with sign change
    % halving interval each loop
    while (R - L)/2 > tol
        M = (L + R)/2;          % midpoint
        if f(L)*f(M) < 0
            R = M;              % root in left half
        else
            L = M;              % root in right hlf
        end
    end
    root = (L + R)/2;           % final midpoint = root
end

function root = newton_method(f, df, x0, tol)
    % fast convergence near root
    % using slope to jump closer
    x = x0;
    while abs(f(x)) > tol
        x = x - f(x)/df(x);     % updating step
    end
    root = x;                   % final x = root
end

function root = secant_method(f, x0, x1, tol)
    % no derivative needed
    % using secant line slope
    while abs(f(x1)) > tol
        x_next = x1 - f(x1)*(x1 - x0)/(f(x1) - f(x0)); % updating step
        x0 = x1;              % shift points
        x1 = x_next;
    end
    root = x1;                % final x = root
end











%% part 2
f = @(x) (x.^3)/100 - (x.^2)/8 + 2*x + 6*sin(x/2 + 6) - 0.7 - exp(x/6);
df = @(x) 3*(x.^2)/100 - x/4 + 2 + 3*cos(x/2 + 6) - exp(x/6)/6;

% where f(l) and f(r) have opposite signs
L = 0; R = 1;

% checking signs to make sure a root is inside
fprintf('f(L) = %.4f, f(R) = %.4f\n', f(L), f(R));

figure;
fplot(f, [-15 40]);     
yline(0,'--');           % horizontal line at y=0
grid on;
title('checking where roots are');
xlabel('x'); ylabel('f(x)');

% midpoint for newton start
M = (L + R)/2;

% running all mthds using new generic solvers
[root_bisect, itb, flagb] = bisection_solver(f, L, R, 1e-6, 1e-6, 200);
[root_newton, itn, flagn] = newton_solver(@(x) deal(f(x), df(x)), M, 1e-6, 1e-6, 100);
[root_secant, its, flags] = secant_solver(f, L, R, 1e-6, 1e-6, 100);

% printing results
fprintf('bisection root: %.6f  iters:%d flag:%d\n', root_bisect, itb, flagb);
fprintf('newton root: %.6f    iters:%d flag:%d\n', root_newton, itn, flagn);
fprintf('secant root: %.6f    iters:%d flag:%d\n', root_secant, its, flags);

function [root, it, flag] = bisection_solver(f, L, R, Athresh, Bthresh, maxit)
    % part 2: bisection safeguard (make sure root is bracketed)
    if f(L)*f(R) > 0
        root = NaN; it = 0; flag = -1; % -1 = bad bracket
        return
    end
    % halving interval each loop
    it = 0; flag = 0;
    while (R - L)/2 > Athresh && it < maxit
        M = (L + R)/2;          
        if abs(f(M)) < Bthresh
            root = M; flag = 1; return   % success by |f(x)| small
        end
        if f(L)*f(M) < 0
            R = M;              % root in left half
        else
            L = M;              % root in right hlf
        end
        it = it + 1;
    end
    root = (L + R)/2; flag = 1; % final midpoint = root
end

function [root, it, flag] = newton_solver(fun_both, x0, Athresh, Bthresh, maxit)
    % fast convergence near root
    % using slope to jump closer
    root = x0; it = 0; flag = 0;
    [fx, dfx] = fun_both(root);
    while abs(fx) > Bthresh && it < maxit
        if dfx == 0
            flag = -2; return  % -2 = zero derivative
        end
        x_new = root - fx/dfx; % updating step
        if abs(x_new - root) < Athresh
            root = x_new; flag = 1; return
        end
        root = x_new;
        [fx, dfx] = fun_both(root);
        it = it + 1;
    end
    flag = 1;
end

function [root, it, flag] = secant_solver(f, x0, x1, Athresh, Bthresh, maxit)
    % no derivative needed
    % using secant line slope
    it = 0; flag = 0;
    f0 = f(x0); f1 = f(x1);
    while abs(f1) > Bthresh && it < maxit
        denom = (f1 - f0);
        if denom == 0
            root = x1; flag = -2; return % -2 = zero denominator
        end
        x_next = x1 - f1*(x1 - x0)/denom; % updating step
        if abs(x_next - x1) < Athresh
            root = x_next; flag = 1; return
        end
        x0 = x1; f0 = f1;
        x1 = x_next; f1 = f(x1);
        it = it + 1;
    end
    root = x1; flag = 1;
end
% output explanation:
% prints f(l) and f(r) to check signs (makes sure root is in interval)
% shows root, number of iterations, and exit flag for each method
% bisection takes more steps since it halves the interval
% newton is fastest if start guess is close, secant is almost as fast
% flags: 1=success, -1=bad bracket, -2=zero division, 0=max iters w/o
% convergence

% part 2 additions:
% athresh = step size cutoff |Δx|, stops when updates are very small
% bthresh = function cutoff |f(x)|, stops when function value is near zero
% maxit   = iteration limit to avoid infinite loops
% solvers now return root, iterations, and flag
% added check for sign change in bisection
% added divide-by-zero checks for newton/secant
% added stop checks for small step size and small f(x)
