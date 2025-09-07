%% part 2 and 3 with step 1-4
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
% step 1: saving iterates in glist1/glist2/glist3
% step 2: roots below are the reference used for errors right after
[root_bisect, itb, flagb, glist1] = bisection_solver(f, L, R, 1e-6, 1e-6, 200);
[root_newton, itn, flagn, glist2] = newton_solver(@(x) deal(f(x), df(x)), M, 1e-6, 1e-6, 100);
[root_secant, its, flags, glist3] = secant_solver(f, L, R, 1e-6, 1e-6, 100);

% printing results
fprintf('bisection root: %.6f  iters:%d flag:%d\n', root_bisect, itb, flagb);
fprintf('newton root: %.6f    iters:%d flag:%d\n', root_newton, itn, flagn);
fprintf('secant root: %.6f    iters:%d flag:%d\n', root_secant, its, flags);

% step 3: random trials to collect (xn, xn+1) pairs
rng(0);                  % fixed seed for repeatability
ntrials = 40;            % how many random runs per method

% containers for all pairs
bef_b = []; aft_b = [];
bef_n = []; aft_n = [];
bef_s = []; aft_s = [];

% bisection trials (random brackets near the same root)
for k = 1:ntrials
    a = -0.2 + 1.4*rand;           % random left in [-0.2, 1.2]
    b = a + 0.6 + 0.8*rand;        % random right at least 0.6 away
    if f(a)*f(b) > 0, continue, end

    [rk, itk, flagk, glk] = bisection_solver(f, a, b, 1e-6, 1e-6, 200);
    if flagk ~= 1 || numel(glk) < 2, continue, end

    % step 1: xn sequence lives in glk
    % step 2: error uses this run's root as reference
    ek = abs(glk - rk);

    % build (eps_n, eps_{n+1})
    bef_b = [bef_b, ek(1:end-1)];
    aft_b = [aft_b, ek(2:end)];
end

% newton trials (random single starts near the root)
for k = 1:ntrials
    x0 = -0.5 + 5.0*rand;          % random start in [-0.5, 1.5]
    [rk, itk, flagk, glk] = newton_solver(@(x) deal(f(x), df(x)), x0, 1e-6, 1e-6, 100);
    if flagk ~= 1 || numel(glk) < 2, continue, end

    ek = abs(glk - rk);
    bef_n = [bef_n, ek(1:end-1)];
    aft_n = [aft_n, ek(2:end)];
end

% secant trials (random pairs near the root)
for k = 1:ntrials
    x0 = -0.7 + 1*rand;
    x1 = -0.7 + 5.0*rand;
    if x0 == x1, x1 = x1 + 1e-3; end

    [rk, itk, flagk, glk] = secant_solver(f, x0, x1, 1e-6, 1e-6, 100);
    if flagk ~= 1 || numel(glk) < 2, continue, end

    ek = abs(glk - rk);
    bef_s = [bef_s, ek(1:end-1)];
    aft_s = [aft_s, ek(2:end)];
end

% quick check that step 3 collected data (and step 1/2/4 are in play)
fprintf('pairs collected -> bisection:%d  newton:%d  secant:%d\n', ...
    numel(bef_b), numel(bef_n), numel(bef_s));
if isempty(bef_b) || isempty(bef_n) || isempty(bef_s)
    warning('some method did not collect data correctly');
else
    disp('all methods produced data -> step 3 good');
end

%example for how to filter the error data
%currently have error_list0, error_list1, index_list
%data points to be used in the regression
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

% combined plot to see all three on one figure
figure; grid on
loglog(bef_b, aft_b, 'ro', 'markerfacecolor','r', 'markersize',3)   % bisection
hold on;
loglog(bef_n, aft_n, 'go', 'markerfacecolor','g', 'markersize',3)   % newton
loglog(bef_s, aft_s, 'bo', 'markerfacecolor','b', 'markersize',3)   % secant
xlabel('\epsilon_n'); ylabel('\epsilon_{n+1}');
title('error map from random trials')
legend('bisection','newton','secant','location','best')
hold off

figure; grid on
loglog(bregx, bregy, 'ro', 'markerfacecolor','r', 'markersize',3)
hold on
loglog(nregx, nregy, 'ro', 'markerfacecolor','g', 'markersize',3)
loglog(sregx, sregy, 'ro', 'markerfacecolor','b', 'markersize',3)
[p_b, k_b] = generate_error_fit(bregx, bregy)
[p_n, k_n] = generate_error_fit(nregx, nregy)
[p_s, k_s] = generate_error_fit(sregx, sregy)
[dfdx,d2fdx2] = approximate_derivative(f,root_newton);
newtexp = abs((1/2)*(d2fdx2/dfdx))



figure
ms = glist1;                                  % step 1: xn sequence (bisection midpoints)
ms_err = abs(ms - root_bisect);               % step 2: error vs chosen reference root
bef = ms_err(1:end-1);                        % step 4: eps_n
aft = ms_err(2:end);                          % step 4: eps_{n+1}
loglog(bef, aft, "r.")                        % step 4: log-log plot of (eps_n, eps_{n+1})
loglog(bef, aft, 'ro', 'markerfacecolor','r','markersize',1);  % step 4: same points styled

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
        glist(end+1) = root;                % step 1: save iterate (kept as is)
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


% part 3 changes compared to part 2
% step 1
% kept saving iterates in glist for each method to capture xn sequence
% no refactor of solvers just made sure glist is recorded every iteration
%
% step 2
% still use the computed root from each run as the error reference
% errors are built as abs(xn - root) for each run
%
% step 3
% added random trials using rand to generate many starts and brackets near the root
% bisection uses random a and b with opposite signs
% newton uses random x0
% secant uses random x0 and x1
% only keep successful runs where flag is 1 and at least two iterates exist
% collect eps_n and eps_n+1 pairs by concatenating across trials
% added rng(0) for repeatability
%
% step 4
% built eps_n and eps_n+1 from the saved glist values
% made loglog plots of eps_n and eps_n+1 to visualize convergence order
% kept plotting style simple to match earlier code


% part 3 additions summary
% added rng(0) so results are repeatable when using random brackets/starts
% added ntrials loop to run many random bisection newton and secant trials
% for each trial collected the iterate list glk and built errors |xn - root|
% concatenated errors across trials to build large arrays of eps_n eps_n+1
% printed how many pairs each method produced as a quick sanity check
% added a combined plot with loglog axes to compare convergence visually
% red = bisection (linear slope ~1)
% green = newton (quadratic slope ~2)
% blue = secant (superlinear slope ~1.6)
% seeing a diagonal cluster for each color means step 3 is working
% bisection gives the most points since it takes more iterations
% newton converges fastest so has fewer points but a steeper slope
% secant sits between bisection and newton in both point count and slope


% note on random trial ranges
% a and b are picked near the root so bisection brackets the same root every time
% interval length kept between 0.6 and 1.4 so it is not too small or too wide
% newton x0 is picked near the root so we see convergence but avoid wrong root
% secant uses two starts in same range and forces x0 ≠ x1 to avoid division by zero
% these ranges give enough variation to get good error pairs without solver failure
