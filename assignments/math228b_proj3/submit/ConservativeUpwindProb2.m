% Generalization of the linear upwind method presented,
% numerically conservative, solving the Burger's Equation
% in project 3, problem 2
% Math 228B, Hongli Zhao (Bob Johnson) - honglizhaobob@berkeley.edu
% exactly same except for the soving U part
% =========== paramters

prob = 2; % 1 - a1, rarefaction solution; 2 - a2, unique shock sol
if prob == 1
    % prob a1
    uL = 0;
    uR = 1;
else
    % prob a2
    uL = 1;
    uR = 0;
end

NX = [60,120,240];

% store error
Trial_errors = zeros(1,3);
Weak_errors = Trial_errors;
% discretizations, change if needed for CFL

for trial = 1:3
    % space discretization
    Nx = NX(trial);
    dx = 6/Nx;
    X = linspace(-3.0,3.0,Nx);

    % time discretization
    FinalTime = 0.8; 
    Nt = Nx;
    dt = FinalTime/Nt;  % change for CFL if needed
    T = linspace(0,FinalTime,Nt);

    % the actual exact weak solution
    Uweak = exactSolution(prob,X,T);
    
    % store solution data
    U = zeros(Nx,Nt); % numerical solution U(xj,tn)

    % boundary data
    U(1,:) = uL; U(Nx,:) = uR;

    % initial data
    U(1:Nx/2,1) = uL; U((Nx/2)+1:Nx,1) = uR;

    % solve directly using upwind
    for nt = 1:Nt-1 % for time
        for nx = 2:Nx-1 % for space
            % fix t solve x first
            U(nx,nt+1) = U(nx,nt) - (dt/dx)*...
                (Flux(U(nx,nt)) - Flux(U(nx-1,nt)));
        end
    end

    % ============================================================
    %   error analysis 1, with base solution
    % ============================================================
    % get "exact" solution
    [UExact,Xexact,Texact] = upwindExact(uL,uR);

    % error matrix
    % use interp2 to fill in the gaps in numerical sols
    [X_grid,T_grid] = meshgrid(X,T);
    [Xexact_grid,Texact_grid] = meshgrid(Xexact,Texact);
    %U_exact_interp = interp2(Xexact_grid,Texact_grid,...
    %   UExact,X_grid,T_grid,'makima');

    % fill in the gaps between two discretizations
    U_interp = interp2(X_grid,T_grid,...
       U,Xexact_grid,Texact_grid,'makima');

    % compute inf error
    disp('====================');
    inf_error = max(max(abs(U_interp - UExact)));
    
    % store error
    Trial_errors(trial) = inf_error;
    
    disp(['Using grid size N = ', num2str(Nx)]);
    disp(['Error measured in inf norm: ', num2str(inf_error)]);
    disp('====================');

    
    
    % can look at what happened over time
    %for row = 1:Nt
    %    Uplot = U(:,row);
    %    figure(1)
    %    plot(X,Uplot,'r','LineWidth',1.5)
    %    title("Numerical Wave Behavior Overtime")
    %end
    Uplot = U(:,Nt);
   
    grid on
    figure(1)
    plot(X,Uplot,'b','LineWidth',1.5)
    title(['Numerical Wave Behavior Overtime',' for N = ', num2str(Nx)])
    xlabel('x')
    ylabel('u(x,FinalTime)')
    grid on

    % ============================================================
    %   error analysis 2, convergence to the actual solution
    %   compare numerical solution with the exact solution
    % ============================================================

    % Error matrix 2 for error from the actual weak solutions
    inf_error_weak = max(max(abs(U-Uweak)));
    
    % store error
    Weak_errors(trial) = inf_error_weak;
    
    disp(['Using grid size N = ', num2str(Nx)]);
    disp(['Error from Weak Solution measured in inf norm: ', ...
        num2str(inf_error_weak)]);
    disp('====================');


    % ============================================================
    %   error analysis 4 (2(d)), at FinalTime=0.8
    %   point wise difference between numerical and exact solution
    % ============================================================
    Uweak_final = Uweak(:,Nt); U_final = U(:,Nt);
    diff = -(Uweak_final - U_final); % pointwise difference
    absdiff = abs(Uweak_final - U_final); % absolute difference
    figure(2)
    subplot(2,1,1)
    plot(X,diff,'b','LineWidth',1.5)
    title('point-wise difference using N = 240, T = 0.8');
    xlabel('x'); ylabel('diff');
    legend('difference behavior')
    grid on
    
    subplot(2,1,2)
    plot(X,absdiff,'r','LineWidth',1.4)
    title('point-wise abs difference using N = 240, T = 0.8');
    xlabel('x'); ylabel('abs diff');
    legend('abs difference behavior')
    grid on
end

% Error plotting...
% log-log plot of base error

disp('...plotting errors')
figure(3)
subplot(1,1,1)
loglog(NX(1:2),6./NX(1:2),NX,Trial_errors,'LineWidth',2.5)
xlabel('num steps N')
ylabel('error')
title("Error Behavior as a Function of N")
grid on
disp('...   done')

% ============================================================
%                       error analysis 3
%   compare numerical solution with the exact weak solution
% ============================================================
disp('...now displaying plot comparison at T=0.8')


subplot(1,1,1)
plotObj = plot(X,Uplot,X,Uweak(:,Nt));
plotObj(1).Color = 'r'; plotObj(1).LineWidth = 1.2;
plotObj(2).Color = 'b'; plotObj(2).LineWidth = 1.2;
title("Spatial Plot Comparison Numerical v Exact")
legend('Numerical','Exact')
xlabel('x')
ylabel('U(x,0.8)')
grid on
disp('...   done')

% =========== utilities
function F = Flux(u)
    % returns numerical flux F(u)
    F = (1/2)*u^2;
end
function [U_exact,X,T] = upwindExact(uL,uR)
    % generates a base solution for error analysis
    % N = 240 discretization
    
    Nx = 240; % 60, 120, 200
    dx = 6/Nx;
    X = linspace(-3.0,3.0,Nx);
    % time discretization
    FinalTime = 0.8; 
    Nt = Nx;
    dt = FinalTime/Nt;  % change for CFL if needed
    T = linspace(0,FinalTime,Nt);

    % store solution data
    U_exact = zeros(Nx,Nt); % numerical solution U(xj,tn)

    % boundary data
    U_exact(1,:) = uL; U_exact(Nx,:) = uR;

    % initial data
    U_exact(1:Nx/2,1) = uL; U_exact((Nx/2)+1:Nx,1) = uR;

    % solve directly using upwind
    for nt = 1:Nt-1 % for time
        for nx = 2:Nx-1 % for space
            
           % fix t solve x first
            U_exact(nx,nt+1) = U_exact(nx,nt) - (dt/dx)*...
            (Flux(U_exact(nx,nt)) - Flux(U_exact(nx-1,nt)));
        end
    end
end


function [U_exact] = exactSolution(prob,X,T)
    % generates exact solutions given the grid points
    
    % Input:
    %       prob - which problem the exact solution is re-
    %   presenting. Takes values of 1 and 2
    %       1, uL < uR, use rarefaction wave solution
    %       2, uL > uR, unique jump discontinuity solution
    %       
    %       X - discretized spatial grid
    %       T - discretized temporal grid
    Nx = length(X); Nt = length(T); 
    U_exact = zeros(Nx,Nt);
    if prob == 1
        uL = 0;
        uR = 1;
    else
        uL = 1;
        uR = 0;
    end
    % start filling in solutions
    if prob == 1
        % prob 1, rarefaction wave sol
        for nt = 1:Nt
            tval = T(nt);
            for nx = 1:Nx
                xval = X(nx);
                if xval < 0
                    U_exact(nx,nt) = 0;
                elseif xval >= 0 && xval < tval
                    U_exact(nx,nt) = xval/tval;
                else
                    U_exact(nx,nt) = 1;
                end
            end 
        end
    else
        % prob 2, unique shock wave sol
        % Rankine Hugoniot speed for shock wave
        s = (uL+uR)/2;
        for nt = 1:Nt
            tval = T(nt);
            for nx = 1:Nx
                xval = X(nx);
                if xval <= s * tval
                    U_exact(nx,nt) = 1;
                else
                    U_exact(nx,nt) = 0;
                end
            end 
        end 
    end
end
