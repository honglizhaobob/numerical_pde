% Solve the 1D heat equation u_t = \kappa u_{xx} 
% using grid size "m" and timestep multiplier "kmul".
% Integrates until final time "T" and plots each solution.
% 1> FTCS scheme (Forward-Time Centered-Space scheme)
% 2> Crank-Nicolson scheme
% 3> Method of Lines with ODE solvers
%
% UC Berkeley Math 228B, Suncica Canic (canics@berkeley.edu)

m = 100;
T = 0.2;
kmul = 0.5;
nsteps=100;

%You can switch between 3 different approaches here.
method =  'odesolver';  %'ftcs';   %'cn'

switch method
    case 'ftcs'
          heat_ftcs(m, T, kmul);
    case 'cn'
          heat_cn(m, T, kmul);
    case 'odesolver'
          heat_odesolver(m, T);
end

% FTCS scheme (Forward-Time Centered-Space scheme)
function heat_ftcs(m, T, kmul)
    h = 1.0 / (m+1); % m grid size
    x = h * (0:m+1);
    k = kmul*h^2;
    N = ceil(T/k);
    
    u = exp(-(x - 0.25).^2 / 0.1^2) + 0.1 * sin(10 * 2 * pi * x);
    u(1) = 0; u(m+2) = 0; % Dirichlet Boundary Condition u(0) = u(1) = 0
    
    clf(); axis([0, 1, -0.1, 1.1]);  plot(x,u);  xlabel("x"); ylabel("u"); hold on;
    for n = 1:N
        u(2:m+1) = u(2:m+1) + k/h^2 * (u(1:m) - 2*u(2:m+1) + u(3:m+2));
        if mod(n, 10) == 0 
            axis([0, 1, -0.1, 1.1]); plot(x,u);pause(1e-3);
        end
    end
end


% Crank-Nicolson scheme
function heat_cn(m, T, kmul)
    h = 1.0 / (m+1); % m grid size
    x = h * (0:m+1);
    k = kmul*h^2;
    N = ceil(T/k);
    
    u = exp(-(x - 0.25).^2 / 0.1^2) + 0.1 * sin(10 * 2 * pi * x);
    u(1) = 0; u(m+2) = 0;
    
    % Form the matrices in the Crank-Nicolson scheme (Left and right)
    nOnes = ones(m, 1) ;
    A = (diag(2 * nOnes, 0) - diag(nOnes(1:m-1), -1) - diag(nOnes(1:m-1), 1))/h^2;

    LH = diag(nOnes, 0) + A*k/2;
    RH = diag(nOnes, 0) - A*k/2;
    clf(); axis([0, 1, -0.1, 1.1]); plot(x,u); hold on;
    for n = 1:N
        u(2:m+1) = LH \ (RH * u(2:m+1)'); 
        if mod(n, 10) == 0 
            axis([0, 1, -0.1, 1.1]); plot(x,u);pause(1e-3);
        end
   end
end

%Matlab ODE solver, see more references at https://www.mathworks.com/help/matlab/math/choose-an-ode-solver.html
function heat_odesolver(m, T)
    h = 1.0 / (m+1);
    x = h * (0:m+1);
    u0 = exp(-(x - 0.25).^2 / 0.1^2) + 0.1 * sin(10 * 2 * pi * x);
    u0(1) = 0; u0(m+2) = 0;
    u=u0;
    
    [t,u] = ode15s(@(t,u) ([0;u(1:m+1)] - 2*u + [u(2:m+2);0]) / h^2 , [0 T], u0);
    
    clf(); axis([0, 1, -0.1, 1.1]); plot(x,u); %hold on;
    for n = 1:length(u)
        axis([0, 1, -0.1, 1.1]); plot(x,u(n,:)), pause(1e-3); hold on;
    end
end
    