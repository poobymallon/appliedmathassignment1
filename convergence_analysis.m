%INPUTS:
%solver_flag: an integer from 1-4 indicating which solver to use
% 1->Bisection 2-> Newton 3->Secant 4->fzero
%fun: the mathematical function that we are using the
% solver to compute the root of
%x_guess0: the initial guess used to compute x_root
%guess_list1: a list of initial guesses for each trial
%guess_list2: a second list of initial guesses for each trial
% if guess_list2 is not needed, then set to zero in input
%filter_list: a list of constants used to filter the collected data
function convergence_analysis(solver_flag, fun, ...
x_guess0, guess_list1, guess_list2, filter_list)

f = @(x) (x.^3)/100 - (x.^2)/8 + 2*x + 6*sin(x/2 + 6) - 0.7 - exp(x/6);
df = @(x) 3*(x.^2)/100 - x/4 + 2 + 3*cos(x/2 + 6) - exp(x/6)/6;

if solver_flag == 1

% this is basically what guess list 1 and 2 should be
% % bisection trials (random brackets near the same root)
% for k = 1:ntrials
%     a = -0.2 + 1.4*rand;           % random left in [-0.2, 1.2]
%     b = a + 0.6 + 0.8*rand;        % random right at least 0.6 away

    % where f(l) and f(r) have opposite signs
    L = x_guess0(1); R = x_guess0(2);
    [root, itb, flagb, glist] = bisection_solver(f, L, R, 1e-14, 1e-14, 500);
    % printing results
    fprintf('bisection root: %.6f  iters:%d flag:%d\n', root, itb, flagb);
    
    bef_b = []; aft_b = [];
    for k = 1:length(guess_list1)
        a = guess_list1(k);
        b = guess_list2(k);
        if f(a)*f(b) > 0, continue, end
    
        [rk, itk, flagk, glk] = bisection_solver(f, a, b, 1e-14, 1e-14, 500);
        if flagk ~= 1 || numel(glk) < 2, continue, end
    
        % step 1: xn sequence lives in glk
        % step 2: error uses this run's root as reference
        ek = abs(glk - rk);
    
        % build (eps_n, eps_{n+1})
        bef_b = [bef_b, ek(1:end-1)];
        aft_b = [aft_b, ek(2:end)];
    end


    bregx = []; % e_n
    bregy = []; % e_{n+1}
    %iterate through the collected data
    for n=1:length(bef_b)
        %if the error is not too big or too small
        %and it was enough iterations into the trial...
        if bef_b(n)>1e-15 && bef_b(n)<1e-2 && ...
                aft_b(n)>1e-14 && aft_b(n)<1e-2 
            %then add it to the set of points for regression
            bregx(end+1) = bef_b(n);
            bregy(end+1) = aft_b(n);
        end
    end
    % quick check that step 3 collected data (and step 1/2/4 are in play)
    fprintf('pairs collected -> bisection:%d', ...
        numel(bef_b));
    if isempty(bef_b) 
        warning('some method did not collect data correctly');
    else
        disp('all methods produced data -> step 3 good');
    end
    [p_b, k_b] = generate_error_fit(bregx, bregy)
end

if solver_flag == 2
% % this is basically what guess list 1 should be 
% % newton trials (random single starts near the root)
% x = -0.5 + 5.0*rand;           
    [root, itn, flagn, glist] = newton_solver(@(x) deal(f(x), df(x)), x_guess0, 1e-14, 1e-14, 500);
    % printing results
    fprintf('newton root: %.6f    iters:%d flag:%d\n', root, itn, flagn);

    bef_n = []; aft_n = [];
    for k = 1:length(guess_list1)
        % random start in [-0.5, 1.5]
        [rk, itk, flagk, glk] = newton_solver(@(x) deal(f(x), df(x)), guess_list1(k), 1e-14, 1e-14, 500);
        if flagk ~= 1 || numel(glk) < 2, continue, end

        ek = abs(glk - rk);
        bef_n = [bef_n, ek(1:end-1)];
        aft_n = [aft_n, ek(2:end)];
    end
    nregx = []; % e_n
    nregy = []; % e_{n+1}
    %iterate through the collected data
    for n=1:length(bef_n)
        %if the error is not too big or too small
        %and it was enough iterations into the trial...
        if bef_n(n)>1e-15 && bef_n(n)<1e-2 && ...
                aft_n(n)>1e-14 && aft_n(n)<1e-2 
            %then add it to the set of points for regression
            nregx(end+1) = bef_n(n);
            nregy(end+1) = aft_n(n);
        end
    end
    % quick check that step 3 collected data (and step 1/2/4 are in play)
    fprintf('pairs collected -> newton:%d', ...
        numel(bef_n));
    if isempty(bef_n)
        warning('some method did not collect data correctly');
    else
        disp('all methods produced data -> step 3 good');
    end
    [p_n, k_n] = generate_error_fit(nregx, nregy)
    [dfdx,d2fdx2] = approximate_derivative(f,root);
    newtexp = abs((1/2)*(d2fdx2/dfdx))
end

if solver_flag == 3
% % this is basically what guess list 1 and 2 should be 
% % secant trials (random pairs near the root)

    [root, its, flags, glist] = secant_solver(f, x_guess0(1), x_guess0(2), 1e-14, 1e-14, 500);
    fprintf('secant root: %.6f    iters:%d flag:%d\n', root, its, flags);

    % containers for all pairs

    bef_s = []; aft_s = [];    
    for k = 1:length(guess_list1)
        x0 = guess_list1(k);
        x1 = guess_list2(k);
        if x0 == x1, x1 = x1 + 1e-7; end
    
        [rk, itk, flagk, glk] = secant_solver(f, x0, x1, 1e-14, 1e-14, 500);
        if flagk ~= 1 || numel(glk) < 2, continue, end
    
        ek = abs(glk - rk);
        bef_s = [bef_s, ek(1:end-1)];
        aft_s = [aft_s, ek(2:end)];
    end
    sregx = []; % e_n
    sregy = []; % e_{n+1}
    %iterate through the collected data
    for n=1:length(bef_s)
        %if the error is not too big or too small
        %and it was enough iterations into the trial...
        if bef_s(n)>1e-15 && bef_s(n)<1e-2 && ...
                aft_s(n)>1e-14 && aft_s(n)<1e-2 
            %then add it to the set of points for regression
            sregx(end+1) = bef_s(n);
            sregy(end+1) = aft_s(n);
        end
    end
    % quick check that step 3 collected data (and step 1/2/4 are in play)
    fprintf('pairs collected ->  secant:%d\n', ...
        numel(bef_s));
    if isempty(bef_s)
        warning('some method did not collect data correctly');
    else
        disp('all methods produced data -> step 3 good');
    end
    [p_s, k_s] = generate_error_fit(sregx, sregy)
end

if solver_flag == 4
    fzero(f, x_guess0)
end


figure;
fplot(fun, [-15 40]);     
yline(0,'--');           % horizontal line at y=0
grid on;
title('checking where roots are');
xlabel('x'); ylabel('f(x)');

%example for how to filter the error data
%currently have error_list0, error_list1, index_list
%data points to be used in the regression



% combined plot to see all three on one figure
figure; grid on
if solver_flag == 1
    loglog(bef_b, aft_b, 'ro', 'markerfacecolor','r', 'markersize',3)   % bisection
end

if solver_flag == 2
    loglog(bef_n, aft_n, 'go', 'markerfacecolor','g', 'markersize',3)   % newton
end

if solver_flag == 3
    loglog(bef_s, aft_s, 'bo', 'markerfacecolor','b', 'markersize',3)   % secant
end
xlabel('\epsilon_n'); ylabel('\epsilon_{n+1}');
title('error map from random trials')
which = {'bisection','newton','secant'};
legend(which(solver_flag),'location','best')
hold off

figure; grid on
if solver_flag == 1
loglog(bregx, bregy, 'ro', 'markerfacecolor','r', 'markersize',3)
end
if solver_flag == 2
loglog(nregx, nregy, 'ro', 'markerfacecolor','g', 'markersize',3)
end
if solver_flag == 3
loglog(sregx, sregy, 'ro', 'markerfacecolor','b', 'markersize',3)
end






figure
ms = glist;                                   % step 1: xn sequence (bisection midpoints)
ms_err = abs(ms - root);               % step 2: error vs chosen reference root
bef = ms_err(1:end-1);                        % step 4: eps_n
aft = ms_err(2:end);                          % step 4: eps_{n+1}
loglog(bef, aft, "r.")                        % step 4: log-log plot of (eps_n, eps_{n+1})
loglog(bef, aft, 'ro', 'markerfacecolor','r','markersize',1);  % step 4: same points styled
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


function [p,k] = generate_error_fit(x_regression,y_regression)
    %generate Y, X1, and X2
    %note that I use the transpose operator (’)
    %to convert the result from a row vector to a column
    %If you are copy-pasting, the ’ character may not work correctly
    Y = log(y_regression)';
    X1 = log(x_regression)';
    X2 = ones(length(X1),1);
    %run the regression
    coeff_vec = regress(Y,[X1,X2]);
    %pull out the coefficients from the fit
    p = coeff_vec(1);
    k = exp(coeff_vec(2));
end


function [dfdx,d2fdx2] = approximate_derivative(fun,x)
%set the step size to be tiny
    delta_x = 1e-6;
%compute the function at different points near x
    f_left = fun(x-delta_x);
    f_0 = fun(x);
    f_right = fun(x+delta_x);
%approximate the first derivative
    dfdx = (f_right-f_left)/(2*delta_x);
%approximate the second derivative
    d2fdx2 = (f_right-2*f_0+f_left)/(delta_x^2);
end
