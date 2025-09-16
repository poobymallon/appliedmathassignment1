function the_egg()
    %set the oval hyper-parameters
    egg_params = struct();
    egg_params.a = 1.5; egg_params.b = 1; egg_params.c = .15;
    %specify the position and orientation of the egg
    x0 = 5; y0 = 5; theta = pi/6;
    %set up the axis
    hold on; axis equal; axis square
    axis([0,10,0,10])
    %plot the origin of the egg frame
    plot(x0,y0,'ro','markerfacecolor','r');
    %compute the perimeter of the egg
    [V_list, ~] = egg_func(linspace(0,1,100),x0,y0,theta,egg_params);
    %plot the perimeter of the egg
    plot(V_list(1,:),V_list(2,:),'k', 'LineWidth', 2);
    %compute a single point along the egg (s=.8)
    %as well as the tangent vector at that point
    [V_single, G_single] = egg_func(.8,x0,y0,theta,egg_params);
    %plot this single point on the egg
    plot(V_single(1),V_single(2),'ro','markerfacecolor','r');
    %plot this tangent vector on the egg
    vector_scaling = .1;
    tan_vec_x = [V_single(1),V_single(1)+vector_scaling*G_single(1)];
    tan_vec_y = [V_single(2),V_single(2)+vector_scaling*G_single(2)];
    plot(tan_vec_x,tan_vec_y,'g')
    [xmin, xmax, ymin, ymax] = compute_bounding_box(x0, y0, theta, egg_params);

    plot([xmin, xmax, xmax, xmin, xmin], [ymin, ymin, ymax, ymax, ymin], 'k--')
    
    traj_params.g = 9.81;
    traj_params.vy = 0;
    traj_params.h0 = 15;
    traj_params.vx = 40;
    traj_params.x0 = 0;
    traj_params.omega = 50;
    traj_params.w0 = 0;
    y_ground = 0;
    x_wall = 20;

    figure;
    hold on;
    [t_ground,t_wall] = collision_func(@egg_trajectory, traj_params, egg_params, y_ground, x_wall);
    t_stop = min(t_ground, t_wall);
    [x0,y0,theta] = egg_trajectory(t_stop, traj_params);
    [V_list, ~] = egg_func(linspace(0,1,100),x0,y0,theta,egg_params);
    %plot the perimeter of the egg
    plot(V_list(1,:),V_list(2,:),'k', 'LineWidth', 2);
    plot([0, x_wall], [y_ground, y_ground], 'k')
    plot([x_wall, x_wall], [0, 10], 'r')
    axis equal;
    
    % animation_example(traj_params,egg_params, x_wall, y_ground, t_stop)
    video_example(traj_params, egg_params, x_wall, y_ground, t_stop)
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

function [V, G] = egg_func(s,x0,y0,theta,egg_params)
    %unpack the struct
    a=egg_params.a;
    b=egg_params.b;
    c=egg_params.c;
    %compute x (without rotation or translation)
    x = a*cos(2*pi*s);
    %useful intermediate variable
    f = exp(-c*x/2);
    %compute y (without rotation or translation)
    y = b*sin(2*pi*s).*f;
    %compute the derivatives of x and y (without rotation or translation)
    dx = -2*pi*a*sin(2*pi*s);
    df = (-c/2)*f.*dx;
    dy = 2*pi*b*cos(2*pi*s).*f + b*sin(2*pi*s).*df;
    %rotation matrix corresponding to theta
    R = [cos(theta),-sin(theta);sin(theta),cos(theta)];
    %compute position and gradient for rotated + translated oval
    V = R*[x;y]+[x0*ones(1,length(theta));y0*ones(1,length(theta))];
    G = R*[dx;dy];
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

function Gx = egg_wrapper_func1_x(s,x0,y0,theta,egg_params)
    [~, G] = egg_func(s,x0,y0,theta,egg_params);
    Gx = G(1);
end

function Gy = egg_wrapper_func1_y(s,x0,y0,theta,egg_params)
    [~, G] = egg_func(s,x0,y0,theta,egg_params);
    Gy = G(2);
end

function [xmin, xmax, ymin, ymax] = compute_bounding_box(x0,y0,theta,egg_params)
    egg_wrapper_func2_x = @(s_in) egg_wrapper_func1_x(s_in,x0,y0,theta,egg_params);
    egg_wrapper_func2_y = @(s_in) egg_wrapper_func1_y(s_in,x0,y0,theta,egg_params);
    s_guesses = 0:.2:1;
    dxtol = 1e-14;
    ytol = 1e-14;
    max_iter = 200;
    dfdxmin = 1e-8;
    x_list = [];
    y_list = [];
    for s_guess = s_guesses
        s_rootx = secant_solver(egg_wrapper_func2_x, s_guess, s_guess+1e-4, dxtol, ytol, max_iter);
        [V,~] = egg_func(s_rootx, x0, y0, theta, egg_params);
        x_list(end+1) = V(1);

        s_rooty = secant_solver(egg_wrapper_func2_y, s_guess, s_guess+1e-4, dxtol, ytol, max_iter);
        [V,~] = egg_func(s_rooty, x0, y0, theta, egg_params);
        y_list(end+1) = V(2);
    end
    xmin = min(x_list);
    xmax = max(x_list);
    ymin = min(y_list);
    ymax = max(y_list);
    % xrange = [xmin, xmax]; 
    % yrange = [ymin, ymax];
end

function [x,y,theta] = egg_trajectory(t, traj_params)
    g = traj_params.g;
    vy = traj_params.vy;
    h0 = traj_params.h0;
    vx = traj_params.vx;
    x0 = traj_params.x0;
    omega = traj_params.omega;
    w0 = traj_params.w0;
    x = vx*t + x0;
    y = -.5*g*t.^2 + vy*t + h0;
    theta = omega*t+w0;
end

function dx = xdist(t, traj_fun, traj_params, egg_params, x_wall)
    [x0, y0, theta] = traj_fun(t, traj_params);
    [~, xmax, ~, ~] = compute_bounding_box(x0,y0,theta,egg_params);
    dx = x_wall - xmax;
end

function dy = ydist(t, traj_fun, traj_params, egg_params, y_ground)
    [x0, y0, theta] = traj_fun(t, traj_params);
    [~, ~, ymin, ~] = compute_bounding_box(x0,y0,theta,egg_params);
    dy = ymin-y_ground;
end

%Function that computes the collision time for a thrown egg
function [t_ground,t_wall] = collision_func(traj_fun, traj_params, egg_params, y_ground, x_wall)
    xdist_wrap = @(t_in) xdist(t_in, traj_fun, traj_params, egg_params, x_wall);
    ydist_wrap = @(t_in) ydist(t_in, traj_fun, traj_params, egg_params, y_ground);
    t_guess = 5;
    t_wall = secant_solver(xdist_wrap, t_guess, t_guess+1e-4, 1e-14, 1e-14, 200);
    t_ground = secant_solver(ydist_wrap, t_guess, t_guess+1e-4, 1e-14, 1e-14, 200);
end

function animation_example(traj_params, egg_params, x_wall, y_ground, t_stop)
    %set up the plotting axis
    figure;
    axis([0, x_wall+1, y_ground, 20])
    hold on; axis equal;
    %initialize the plot
    % Draw static wall and ground
    plot([0, x_wall], [y_ground, y_ground], 'k');  % ground
    plot([x_wall, x_wall], [y_ground, 20], 'r');   % wall
    
    % Initialize the egg plot with empty data
    start_plot = plot(nan, nan, 'k', 'LineWidth', 2);
    %iterate through time
    for t=0:.001:t_stop
        [x0,y0,theta] = egg_trajectory(t, traj_params);
        [V_list, ~] = egg_func(linspace(0,1,100),x0,y0,theta,egg_params);
        % %plot the perimeter of the egg
        % plot(V_list(1,:),V_list(2,:),'k', 'LineWidth', 2);
        % plot([0, x_wall], [y_ground, y_ground], 'k')
        % plot([x_wall, x_wall], [0, 10], 'r')
        %update the coordinates of the square plot
        set(start_plot,'xdata',V_list(1,:),'ydata',V_list(2,:));
        %update the actual plotting window
        drawnow;
    end
end

function video_example(traj_params, egg_params, x_wall, y_ground, t_stop)
    %define location and filename where video will be stored
    %written a bit weird to make it fit when viewed in assignment
    mypath = 'C:\orionmath\appliedmathassignment1\';
    fname='egg_animation.avi';
    input_fname = [mypath,fname];
    %create a videowriter, which will write frames to the animation file
    writerObj = VideoWriter(input_fname);
    %must call open before writing any frames
    open(writerObj);

    %set up the plotting axis
    clf
    fig1 = figure(3);
    axis([0, x_wall+1, y_ground, 20])
    hold on; axis equal;
    %initialize the plot
    % Draw static wall and ground
    plot([0, x_wall], [y_ground, y_ground], 'k');  % ground
    plot([x_wall, x_wall], [y_ground, 20], 'r');   % wall
    
    % Initialize the egg plot with empty data
    start_plot = plot(nan, nan, 'k', 'LineWidth', 2);
    %iterate through time
    for t=0:.001:t_stop
        [x0,y0,theta] = egg_trajectory(t, traj_params);
        [V_list, ~] = egg_func(linspace(0,1,100),x0,y0,theta,egg_params);
        % %plot the perimeter of the egg
        % plot(V_list(1,:),V_list(2,:),'k', 'LineWidth', 2);
        % plot([0, x_wall], [y_ground, y_ground], 'k')
        % plot([x_wall, x_wall], [0, 10], 'r')
        %update the coordinates of the square plot
        set(start_plot,'xdata',V_list(1,:),'ydata',V_list(2,:));
        %update the actual plotting window
        drawnow;

        %capture a frame (what is currently plotted)
        current_frame = getframe(fig1);
        %write the frame to the video
        writeVideo(writerObj,current_frame);
    end
    for i = 1:20
        [x0,y0,theta] = egg_trajectory(t_stop, traj_params);
        [V_list, ~] = egg_func(linspace(0,1,100),x0,y0,theta,egg_params);
        % %plot the perimeter of the egg
        % plot(V_list(1,:),V_list(2,:),'k', 'LineWidth', 2);
        % plot([0, x_wall], [y_ground, y_ground], 'k')
        % plot([x_wall, x_wall], [0, 10], 'r')
        %update the coordinates of the square plot
        set(start_plot,'xdata',V_list(1,:),'ydata',V_list(2,:));
        %update the actual plotting window
        drawnow;

        %capture a frame (what is currently plotted)
        current_frame = getframe(fig1);
        %write the frame to the video
        writeVideo(writerObj,current_frame);
    end
    %must call close after all frames are written
    close(writerObj);


end