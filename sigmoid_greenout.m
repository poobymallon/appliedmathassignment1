function sigmoid_greenout()
    x_in = linspace(0,50,200);
    [fvals, dfdxvals] = test_function03(x_in);

    dxtol = 1e-14;
    ytol = 1e-14;
    dfdxmin = 1e-10;
    max_iter = 100;

    x_root = bisection_solver(@test_function03, 0, 50, 1e-14, 1e-14, 200);
    % x_root = newton_solver(@test_function03, 27, 1e-14, 1e-14, 200);
    % x_root = secant_solver(@test_function03, 26, 27, 1e-14, 1e-14, 200);
    % x_root = fzero(@test_function03, 27);
    
    success_list = [];
    fail_list = [];

    for n = 1:length(x_in)
        x_guess = x_in(n);
        root_attempt = bisection_solver(@test_function03, x_guess-25, x_guess+25, 1e-14, 1e-14, 200);
        % root_attempt = newton_solver(@test_function03, x_guess, 1e-14, 1e-14, 200);
        % root_attempt = secant_solver(@test_function03, x_guess-1, x_guess, 1e-14, 1e-14, 200);
        % root_attempt = fzero(@test_function03, x_guess);

        if abs(x_root - root_attempt) < .1
            success_list(end+1) = x_guess;
        else
            fail_list(end+1) = x_guess;
        end
    end
    figure
    hold on

    [f_success , ~] = test_function03(success_list);
    [f_fail , ~] = test_function03(fail_list);
    plot(x_in, fvals, 'k')
    plot(success_list, f_success, 'g.')
    plot(fail_list, f_fail, 'r.')
end

%Example sigmoid function
function [f_val,dfdx] = test_function03(x)
% global input_list;
% input_list(:,end+1) = x;
a = 27.3; b = 2; c = 8.3; d = -3;
H = exp((x-a)/b);
dH = H/b;
L = 1+H;
dL = dH;
f_val = c*H./L+d;
dfdx = c*(L.*dH-H.*dL)./(L.^2);
end

% function x_root = orion_newton(fun, x_guess, dxtol, ytol, max_iter, dfdxmin)
%     delta_x = 2*dxtol;
%     [fval,dfdx] = fun(x_guess);
%     count = 0;
% 
%     while count<max_iter && abs(delta_x) > dxtol && abs(fval) > ytol && abs(dfdx) > dfdxmin
%         count = count+1;
%         delta_x = -fval/dfdx;
%         x_guess = x_guess + delta_x;
%     end
%     x_root = x_guess;
% end

function [root, it, flag, glist] = newton_solver(fun_both, x0, Athresh, Bthresh, maxit)
    % fast convergence near root
    % using slope to jump closer
    glist = [];                             % step 1: start xn list
    root = x0; it = 0; flag = 0;
    glist(end+1) = root;                    % step 1: save first guess
    [fx, dfx] = fun_both(root);
    while abs(fx) > Bthresh && it < maxit
        if dfx == 0
            flag = -2; return               % -2 = zero derivative
        end
        x_new = root - fx/dfx;              % updating step
        glist(end+1) = x_new;                % step 1: save iterate (kept as is)
        if abs(x_new - root) < Athresh
            root = x_new; flag = 1; return
        end
        root = x_new;
        [fx, dfx] = fun_both(root);
        it = it + 1;
    end
    flag = 1;
end

function [root, it, flag, glist] = secant_solver(f, x0, x1, Athresh, Bthresh, maxit)
    % no derivative needed
    % using secant line slope
    it = 0; flag = 0;
    f0 = f(x0); f1 = f(x1);
    glist = [];                             % step 1: start xn list
    glist(end+1) = x0;                      % step 1: save first two guesses
    glist(end+1) = x1;                      % step 1: save first two guesses
    while abs(f1) > Bthresh && it < maxit
        denom = (f1 - f0);
        if denom == 0
            root = x1; flag = -2; return    % -2 = zero denominator
        end
        x_next = x1 - f1*(x1 - x0)/denom;   % updating step
        glist(end+1) = x_next;              % step 1: save iterate
        if abs(x_next - x1) < Athresh
            root = x_next; flag = 1; return
        end
        x0 = x1; f0 = f1;
        x1 = x_next; f1 = f(x1);
        it = it + 1;
    end
    root = x1; flag = 1;
end

function [root, it, flag, glist] = bisection_solver(f, L, R, Athresh, Bthresh, maxit)
    % part 2: bisection safeguard (make sure root is bracketed)
    if f(L)*f(R) > 0
        root = NaN; it = 0; flag = -1; % -1 = bad bracket
        return
    end
    glist = [];                             % step 1: start xn list
    % halving interval each loop
    it = 0; flag = 0;
    while (R - L)/2 > Athresh && it < maxit
        M = (L + R)/2;  
        glist(end+1) = M;                   % step 1: save iterate (midpoint)
        if abs(f(M)) < Bthresh
            root = M; flag = 1; return      % success by |f(x)| small
        end
        if f(L)*f(M) < 0
            R = M;                          % root in left half
        else
            L = M;                          % root in right hlf
        end
        it = it + 1;
    end
    root = (L + R)/2; flag = 1;             % final midpoint = root
end