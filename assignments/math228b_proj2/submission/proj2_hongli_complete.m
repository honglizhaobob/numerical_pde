%  Math 228B, Hongli Zhao

%  Completed Project 2, solver for Euler's gas equations  
%  u_t +div(F) = 0 with periodic boundary conditions
%  Divergence calculation using 4th order Padé scheme with filtering
%  Numerical integration using 4th order Runga-kutta scheme

% ============================================================
% function calls, organized by Proj2 spec
clear
% (f) Isentropic Vortex Convection problem 
% and Kelvin-Helmholtz instability

% control sequences:
problem = 0; % 0 - vortex, 1 - kelvin
with_plot = 1; % turn with_plot = 1 to see final sol
see_error = 0; % turn see_error = 1 to see error plot

% load parameters given problem 0, 1
if problem == 0
    % vortex problem
    
    % define data holder for plotting
    % 2 (alpha choice) by 3 (step size choice)
    inf_err_data = zeros(2,3); 
    inf_r_err_data = inf_err_data;
    inf_ru_err_data = inf_err_data;
    inf_rv_err_data = inf_err_data;
    inf_rE_err_data = inf_err_data;
    
    % change N for different results
    Ns = [32,64,128]; % N = 32,64,128
    alphas = [0.48,0.499]; 
    
    for alpha_idx=1:2 % over alpha
        for stpsz_idx=1:3 % over step size
            alpha = alphas(alpha_idx);
            % grid
            T = 5*sqrt(2);
            FinalTime = T;
            h = 10.0 / Ns(stpsz_idx);
            x = h * (0:Ns(stpsz_idx));
            y = h * (0:Ns(stpsz_idx));
            [x,y] = meshgrid(x,y);
            k = 0.3*h;
            adjust_factor = ceil(T/k);
            k = T/adjust_factor;
            % parameters
            gamma = 7/5; % gas constant
            b = 0.5; % vortex strength
            x_c = 5; y_c = 5;
            % initial solution
            d = sqrt((x - x_c).^2 + (y - y_c).^2);
            r = (1 - ((b^2)*(gamma-1)/(8*gamma*pi^2))...
                *exp(1-d.^2)).^(1/(gamma-1));
            p = r.^gamma;
            u = 0.1 - (0.5*b/pi)*exp(0.5*(1-d.^2)).*(y-y_c);
            v = (0.5*b/pi)*exp(0.5*(1-d.^2)).*(x-x_c);
            ru = r.*u;
            rv = r.*v;
            rE = p/(gamma-1) + r.*(u.^2 + v.^2)/2;
            
            % prepare exact solution
            u_inf = 0.1; v_inf = 0;
            r_inf = 1; p_inf = 1;
            x_exact = x_c + FinalTime*u_inf;
            d = sqrt((x-x_exact).^2 + (y-y_c).^2);
            r_exact = (1 -((gamma-1)*(b^2)/(8*gamma*pi^2))...
                *exp(1-d.^2)).^(1/(gamma-1));
            u_exact = u_inf - ((1/2)*b/pi)*exp((1/2)...
                *(1-d.^2)).*(y-y_c);
            v_exact = v_inf + ((1/2)*b/pi)*exp((1/2)...
                *(1-d.^2)).*(x-x_exact);
            p_exact = r_exact.^gamma;
            ru_exact = r_exact.*u_exact;
            rv_exact = r_exact.*v_exact;
            rE_exact = p_exact/(gamma-1)...
                +r_exact.*(u_exact.^2 + v_exact.^2)/2;
                  
            %===================== Numerical sol
            for n = 1:adjust_factor
               [r,ru,rv,rE] = euler_rk4step(r,ru,rv,rE,h,k,alpha);
            end

            % Error matrices
            r_error_mat = r - r_exact;
            ru_error_mat = ru - ru_exact;
            rv_error_mat = rv - rv_exact;
            rE_error_mat = rE - rE_exact;
            % display info, comment out if wanted
            statement1 = ['current alpha: ', num2str(alpha)];
            statement2 = ['current stepsize: '...
                , num2str(Ns(stpsz_idx))];
            
            disp("====================");
            disp(statement1);
            disp(statement2);
            
            % calculate infinity norms for errors
            r_error = max(max(abs(r_error_mat)));
            ru_error = max(max(abs(ru_error_mat)));
            rv_error = max(max(abs(rv_error_mat)));
            rE_error = max(max(abs(rE_error_mat)));
            disp("====================");
            disp(['rho error is: ', num2str(r_error)]);
            disp(['rho u error is: ', num2str(ru_error)]);
            disp(['rho v error is: ', num2str(rv_error)]);
            disp(['rho E error is: ', num2str(rE_error)]);
            
            % overall error
            error_total = r_error+ru_error+rv_error+rE_error;
            disp("====================");
            disp(['total error is: ', num2str(error_total)]);
            disp("====================");
            inf_err_data(alpha_idx,stpsz_idx) = error_total;
            inf_r_err_data(alpha_idx,stpsz_idx) = r_error;
            inf_ru_err_data(alpha_idx,stpsz_idx) = ru_error;
            inf_rv_err_data(alpha_idx,stpsz_idx) = rv_error;
            inf_rE_err_data(alpha_idx,stpsz_idx) = rE_error;
        end
    end
    
    % populated inf_err_data, can plot
    if see_error == 1
        H = 10./Ns;
        error_data = zeros(size(inf_r_err_data));
        error_data(1,:) = inf_r_err_data(1,:) ...
            + inf_ru_err_data(1,:) + ...
            inf_rv_err_data(1,:) + ...
            inf_rE_err_data(1,:);
        
        error_data(2,:) = inf_r_err_data(2,:) ...
            + inf_ru_err_data(2,:) + ...
            inf_rv_err_data(2,:) + ...
            inf_rE_err_data(2,:);
        coe1 = polyfit(log(H),log(error_data(1,:)),1);
        slope1 = coe1(1);
        disp(['using \alpha = ',num2str(alphas(1)),...
            ', error is ',num2str(slope1),' th order'])
        coe2 = polyfit(log(H),log(error_data(2,:)),1);
        slope2 = coe2(1);
        disp(['using \alpha = ',num2str(alphas(2)),...
            ', error is ',num2str(slope2),' th order'])
        %============= granular plotting
        figure(1)
        subplot(1,2,1)
        loglog(H,H.^4,H,error_data(1,:),'LineWidth',2.5)
        
        title(['error behavior, \alpha=',num2str(alphas(1))])
        legend('h^4','agg err')
        subplot(1,2,2)
        loglog(H,H.^4,H,error_data(2,:),'LineWidth',2.5)
        title(['error behavior, \alpha=',num2str(alphas(2))])
        legend('h^4','agg err')
    end
else
    % kelvin problem
    N = 256;
    h = 1.0/N;
    x = (0:h:1.0); y = (0:h:1.0);
    T = 1.0;
    k = 0.3*h;
    adjust_factor = ceil(T/k); % factor to adjust k
    k = T/adjust_factor;
    [x,y] = meshgrid(x,y);
    gamma = 7/5; % gas constant
    % initial conditions
    r = zeros(size(x));
    domain1 = find(abs(y-0.5)<(0.15 + sin(2*pi*x)/200));
    domain2 = find(~(abs(y-0.5)<0.15 + sin(2*pi*x)/200));
    r(domain1) = 2;
    r(domain2) = 1;
    u = r-1;
    v = zeros(size(r));
    p = 3*ones(size(r));
    ru = r.*u;
    rv = r.*v;
    rE = p/(gamma-1)+r.*(u.^2 + v.^2)/2;
end

% time integration
% integrate in time for kelvin
if problem == 1
    alf = 0.480;
    for t = 0:k:T
        [r,ru,rv,rE] = euler_rk4step(r,ru,rv,rE,h,k,alf);
        figure(2)
        contourf(x,y,r,16)
        shading interp
        title('\rho_{FinalTime=1.0}')
        colorbar
    end
end

% plotting, with control 0 and 1
if with_plot ~= 0
    if problem == 0
        % see what happened
        figure(1)
        subplot(2,2,1)
        contourf(x,y,r,16)
        shading interp
        title('\rho')
        colorbar

        subplot(2,2,2)
        contourf(x,y,ru,16)
        shading interp
        title('\rho u')
        colorbar

        subplot(2,2,3)
        contourf(x,y,rv,16)
        shading interp
        title('\rho v')
        colorbar

        subplot(2,2,4)
        contourf(x,y,rE,16)
        shading interp
        title('\rho E')
        colorbar
    else
        figure(2)
        contourf(x,y,r,16)
        shading interp
        title('\rho_{FinalTime=1.0}')
        colorbar
    end
end


% ============================================================
% utilities

function divF = compact_div(Fx,Fy,h)
    % calculates the divergence of a grid function field
    % using the fourth order Padé method with periodic boundary
    % conditions

    % input:
    %   - Fx: (N+1)x(N+1) matrix; the flux function in the x-direction
    %   - Fy: (N+1)x(N+1) matrix; the flux function in the y-direction
    %   - h: grid spacing h

    % output:
    %   - divF: the matrix representing the divergence, for
    % a pair of components e.g. (rx,ry), it is a matrix resulted
    % from adding up two matrices.

    % preallocate matrices to store our divergence in each
    % component. In DX, we compute row-wise (y is fixed) and 
    % similarly for DY.

    % ============================================================
    DX = zeros(size(Fx)); DY = zeros(size(Fy)); 
    n = length(DX)-1;
    % create coefficient matrix A, we solve Af' = (3/h)*f
    A = eye(n,n) * 4;
    for i=1:n-1
       A(i+1,i) = 1; % subdiag
       A(i,i+1) = 1; % superdiag
    end
    A(1,end) = 1; A(end,1) = 1;
    A = sparse(A); % sparsify
    
    % create rhs coefficient matrix B, and apply on data to 
    % create b
    B = zeros(n,n);
    for i=1:n-1
       B(i+1,i) = -1; % subdiag
       B(i,i+1) = 1; % supdiag
    end
    B(1,end) = -1; B(end,1) = 1; % corner
    B = sparse(B); % make sparse
    
    f_data = Fx(1:n+1,1:n)';
    f_data = (3/h) * B * f_data;
    f_prime_data = A\f_data;
    DX(1:n+1,1:n) = f_prime_data';
    
    DX(:,end)=DX(:,1); % periodic


    f_data = Fy(1:n,1:n+1);
    f_data = (3/h) * B * f_data;
    f_prime_data = A\f_data;
    DY(1:n,1:n+1) = f_prime_data;
    
    
    DY(end,:)=DY(1,:); % periodic
    
    % add up to obtain divergence
    divF = DX + DY;
end



function [u] = compact_filter(u, alpha)
    % filters the numerical solution u

    % inputs:
    %   - u: our numerical solution, 4x(n+1)x(n+1) data
    % matrix storing values of r,ru,rv,rE defined on the x-y grid
    % when filtering, do for each component.

    %   - alpha: parameter

    % output:
    %   - u_new: the filtered smooth solution u computed using the
    % scheme, still a 4x(N+1)x(N+1) data matrix.
    % ======================================================
    size_u = size(u); 
    n = size_u(3)-1;
    % compute parameters:
    a = (5/8) + (3*alpha)/4; 
    b = alpha + 1/2; 
    c = alpha/4 - 1/8;

    % tridiag matrix on lhs
    A = eye(n,n);
    for i=1:n-1
        A(i+1,i) = alpha; % subdiag
        A(i,i+1) = alpha; % superdiag
    end
    A(1,end) = alpha; A(end,1) = alpha;
    A = sparse(A); % sparsify
    % prepare right hand side in each loop
    rhs = zeros(n,n+1); % for all components
    for i = 1:4
        % first obtain data, which is (N+1)x(N+1)    
        data = squeeze(u(i,:,:));
        % filter in x
        data(:,n+1) = data(:,1);

        % prepare right hand side
        % first take care of boundaries
        rhs(1,:) = (a*data(:,1)+(c/2)*(data(:,3)+data(:,n-1))...
            + (b/2)*(data(:,2)+data(:,n)))';
        rhs(2,:) = (a*data(:,2)+(c/2)*(data(:,4)+data(:,n))...
            + (b/2)*(data(:,3)+data(:,1)))';
        
        rhs(3:n-1,:) = (a*data(:,3:n-1)+...
            (c/2)*(data(:,5:n+1)+data(:,1:n-3))...
            + (b/2)*(data(:,4:n)+data(:,2:n-2)))';
        rhs(n,:) = (a*data(:,n)+(c/2)*(data(:,2)+data(:,n-2))...
            + (b/2)*(data(:,1)+data(:,n-1)))';
        
        % reassign boundary u(1)=u(N+1)
        data(:,1:n) = (A\rhs)';
        data(:,n+1) = data(:,1);
         
        % filter in y
        data(n+1,:) = data(1,:);
        data(:,n+1) = data(:,1);
        
        % prepare right hand side
        rhs(1,:) = a*data(1,:)+(c/2)*(data(3,:)+data(n-1,:))...
                + (b/2)*(data(2,:)+data(n,:));
        rhs(2,:) = a*data(2,:)+(c/2)*(data(4,:)+data(n,:))...
            + (b/2)*(data(3,:)+data(1,:));
        rhs(3:n-1,:) = a*data(3:n-1,:)+...
            (c/2)*(data(5:n+1,:)+data(1:n-3,:))...
                + (b/2)*(data(4:n,:)+data(2:n-2,:));
        rhs(n,:) = a*data(n,:)+(c/2)*(data(2,:)+data(n-2,:))...
                + (b/2)*(data(1,:)+data(n-1,:));

        % reassign boundary
        data(1:n,:) = A\rhs;
        data(n+1,:) = data(1,:);  
        
        % put back filtered data  
        u(i,:,:) = data;           
    end
    % do four times, each time for one entry, we have 4x(N+1)x(N+1)    
end





function [Frx,Fry,Frux,Fruy,Frvx,Frvy,FrEx,FrEy] = euler_fluxes(r, ...
    ru, rv, rE)
    % inputs: r, ru, rv, rE (vector components of u), each is a scalar
    % (written in general case so easily switched if needed)
    % outputs: flux functions of F, scalar, at (x,y)
    % ======================================================
    gamma = 7/5;
    u = ru./r;
    v = rv./r;
    p = (gamma-1)*(rE - (ru.*u + rv.*v)/2);
    Frx = ru; 
    Fry = rv;
    Frux = ru.*u + p;
    Fruy = ru.*v;
    Frvx = Fruy;
    Frvy = rv.*v + p;
    FrEx = u.*(rE + p);
    FrEy = v.*(rE + p);
end




function [fr,fru,frv,frE] = euler_rhs(r, ru, rv, rE, h)
    % computes the right-hand side of the discretized divergence
    % there is a negative sign

    % inputs:
    %   - r, ru, rv, rE: input (N+1)x(N+1) matrices from which we 
    % would like to calculate the divergence for
    %   - h: step size of your discretization, to be fed in when 
    % approximating divergence

    % outputs: 
    %   - fr, fru, fry, frE: the four divergence components
    % ======================================================

    % Get all 8 fluxes:
    [Frx,Fry,Frux,Fruy,Frvx,Frvy,FrEx,FrEy] = euler_fluxes(r, ...
    ru, rv, rE);

    fr = (-1) * compact_div(Frx,Fry,h);
    fru = (-1) * compact_div(Frux,Fruy,h);
    frv = (-1) * compact_div(Frvx,Frvy,h);
    frE = (-1) * compact_div(FrEx,FrEy,h);

end




function [r,ru,rv,rE] = euler_rk4step(r,ru,rv,rE,h,k,alpha)
    % calculate our final integrated solution
    % using 4th order accurate Runge-Kutta method
    % integrator in time t

    % inputs:
    %   - r, ru, rv, rE: the function values (at each (x,y) on grid)
    %   computed and organized in (N+1)x(N+1) matrices
    %   - h, alpha: parameters to feed into our compact filter
    %   - k: our time step, delta t

    % outputs:
    %   - r, ru, rv, rE: the function values evolved after 1 timestep

    % right hand side rhs is computed using euler_rhs
    % can consider as Ut = F(U), where U = [r,ru,rv,rE], F calcs div
    % ======================================================

    % ======= k1
    [fr1,fru1,frv1,frE1] = euler_rhs(r,ru,rv,rE,h);
    r_k1 = k * fr1; ru_k1 = k * fru1;
    rv_k1 = k * frv1; rE_k1 = k * frE1;

    % ======= k2
    [fr2,fru2,frv2,frE2] = ...
        euler_rhs(r + 0.5 * r_k1,ru + 0.5 * ru_k1, ...
            rv + 0.5 * rv_k1,rE + 0.5 * rE_k1,h);

    r_k2 = k * fr2; ru_k2 = k * fru2;
    rv_k2 = k * frv2; rE_k2 = k * frE2;
    % ======= k3
    [fr3,fru3,frv3,frE3] = ...
        euler_rhs(r + 0.5 * r_k2,ru + 0.5 * ru_k2, ...
            rv + 0.5 * rv_k2,rE + 0.5 * rE_k2,h);
    r_k3 = k * fr3; ru_k3 = k * fru3;
    rv_k3 = k * frv3; rE_k3 = k * frE3;
    % ======= k4
    [fr4,fru4,frv4,frE4] = ...
        euler_rhs(r + r_k3,ru + ru_k3, ...
            rv + rv_k3,rE + rE_k3,h);
    r_k4 = k * fr4; ru_k4 = k * fru4;
    rv_k4 = k * frv4; rE_k4 = k * frE4;

    % weighted average
    r_raw = r + (1/6) * (r_k1 + 2*r_k2 + 2*r_k3 + r_k4);
    ru_raw = ru + (1/6) * (ru_k1 + 2*ru_k2 + 2*ru_k3 + ru_k4);
    rv_raw = rv + (1/6) * (rv_k1 + 2*rv_k2 + 2*rv_k3 + rv_k4);
    rE_raw = rE + (1/6) * (rE_k1 + 2*rE_k2 + 2*rE_k3 + rE_k4);
    
    % in order to use compact_filter, need to assemble into u_raw
    size_u = [4, size(r)];
    u_raw = zeros(size_u);
    
    u_raw(1,:,:) = r_raw;
    u_raw(2,:,:) = ru_raw; 
    u_raw(3,:,:) = rv_raw;
    u_raw(4,:,:) = rE_raw;
    
    % filter our solution in each component
    % returns a 4x(n+1)x(n+1) matrix to dessemble
    u_clean = compact_filter(u_raw,alpha);
    
    % dessemble into our output
    r = reshape(u_clean(1,:,:),size(r));
    ru = reshape(u_clean(2,:,:),size(ru));
    rv = reshape(u_clean(3,:,:),size(rv));
    rE = reshape(u_clean(4,:,:),size(rE));
end