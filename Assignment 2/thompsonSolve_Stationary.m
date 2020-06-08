
function ts = thompsonSolve_Stationary(n,d,epsilon,max_iterations,seed)
    
    rng(seed)
    
    n = n-1;
    n_bar = d*n; % calculate n_bar
    
    % INIT
    % define x0 as n d-length vectors between with values between -1 and 1.
    x0_stuck = [1 repelem(0,d-1)]; % set one of the electrons to a constant point
    x0 = -1 + (1+1).*rand(n,d-1);
    x0_lastcol = repelem(1,n);
    % make sure x0 works with conditions (ie: |x_i| - 1 = 0)
    for j = 1:n
        x0_lastcol(j) = x0_lastcol(j) - sum(x0(j,:).^2);
        x0_lastcol(j) = (randi(0:1)*2-1) * sqrt(abs(x0_lastcol(j))); % randomise positive and negative
    end
    x0 = [x0, transpose(x0_lastcol);x0_stuck]; % add final column to xo
    
    u0 = rand(n,1); % initial values of u distributed uniformly on [0,1]
    
    x = x0;
    u = u0;
    repeat = true; % repeat while statement
    f_vect = []; % vect of f(x) over time
    norm_direction_vect = []; % vect of norm(direction) over time
    
    num = 0;
    while repeat == true
        
        num = num+1;
        
        % f calculations:
        f = func(x(n,:)); % Calculate f
        
        % Numerical calculation of grad of f
        g_f = zeros(n_bar,1);
        bumpg = epsilon; % assign bump value
        l = 0;
        for i1 = 1:n
            for i2 = 1:d
                l = l+1;
                bump_xi = x;
                bump_xi(i1,i2) = bump_xi(i1,i2) + bumpg;
                f_bump_xi = func(bump_xi);
                f_init = func(x);    
                g_f(l) = (f_bump_xi - f_init)/bumpg;
            end
        end
        
        % Numerical calculation of the Hessian of f
        bumpH = epsilon; % assign bump value
        H_f = zeros(n_bar); 
        l = 0;
        for i1 = 1:n
            for i2 = 1:d
                for j1 = 1:n
                    for j2 = 1:d
                        l = l+1;
                        bump_xij = x; bump_xi = x; bump_xj = x;
                        bump_xij(i1,i2) = bump_xij(i1,i2) + bumpH;
                        bump_xij(j1,j2) = bump_xij(j1,j2) + bumpH;
                        bump_xi(i1,i2) = bump_xi(i1,i2) + bumpH;
                        bump_xj(j1,j2) = bump_xj(j1,j2) + bumpH;
                        f_bump_xij = func(bump_xij);
                        f_bump_xi = func(bump_xi);
                        f_bump_xj = func(bump_xj);
                        f_init = func(x);    
                        H_f(l) = (f_bump_xij - f_bump_xi - f_bump_xj + f_init)/bumpH^2;
                    end
                end
            end
        end
        
        f_vect(num) = f; % add present value to f vector
        
        % h calculations:
        hi = zeros(n,1); % Calculate h_i vect
        for i = 1:n
            hi(i) = norm(x(i,:))^2 - 1;
        end
        
        g_h = zeros(n,n_bar); % Calculate grad of h_i
        j = 0;
        for i = 0:n-1
            j = j+1;
            g_h(j,i*d+1:(i+1)*d) = 2*x(i+1,:);
        end
        
        H_h = cell(n,1); % Calculate hessian of h_i
        for i = 1:n
            H_h{i} = zeros(n_bar);
            H_h{i}(((i-1)*d+1):(i*d),((i-1)*d+1:i*d)) = diag(repelem(2,d));
        end
        
        Q = H_f; % Calculate Q = H_f + sum(u_i * H_i)  
        for i = 1:n
            Q = Q + u(i) * H_h{i};
        end
        
        % Calculate direction and new multipliers
        P = [Q transpose(g_h); g_h zeros(n)];
        V = [-g_f;-hi];
        dv_vect = inv(P)*V; % vector of directions and multipliers
        
        % seperate vector into direction vector and multiplier vector
        direction = dv_vect(1:n_bar);
        v = dv_vect(n_bar+1:n_bar+n);
        
        new_direction = vec2mat(direction,d,n);
        new_direction = [new_direction;repelem(0,d)];
        x = x + new_direction; % update x by the direction
        
        u = v; % update multipliers
        
        % set end condition
        disp("norm of direction")
        disp(norm(direction))
        if norm(direction) <= epsilon
            repeat = false;
        end
        
        norm_direction_vect(num) = norm(direction);
        
        % set a max number of repetitions
        if num >= max_iterations
            repeat = false;
        end
    end
    
    % Plot result for the 2 dimensionnal case
    figure(1);
    if d == 2
        hold on
        circle(0,0,1) % draw a circle
        scatter(x(:,1),x(:,2),"red", 'MarkerFaceColor', 'red') % Stable points
        scatter(x0(:,1),x0(:,2),"bla", 'MarkerFaceColor', 'bla') % Initial points
        hold off
    end
    
    % Plot result for the 3 dimensionnal case
    if d == 3
        
            for j = i+1:n
                dij = norm(x(i,:) - x(j,:));
                disp(dij)
            end
        end
        [X,Y,Z] = sphere(20);
        surf(X,Y,Z,'FaceColor', 'none') % plot sphere
        hold on
        scatter3(x(:,1),x(:,2),x(:,3),"red","filled") % Stable points
        %scatter3(x0(:,1),x0(:,2),x0(:,3),"bla","filled") % Initial points
        hold off
    end
    
    
    
    ts = x; % return the position of all the points
    % display the solution: min f(x)
    disp("min f(x) = ");disp(func(x))
    
    % plot f(x) over time
    figure(2);
    plot(f_vect);
    %plot norm(direction) over time
    figure(3);
    plot(norm_direction_vect);
    
end

