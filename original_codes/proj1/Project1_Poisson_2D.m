% Solve 2D Poisson's equation -(u_xx+u_yy) = 1, on trapezoidal domain 
% with Dirichlet (bottom,right) and Neumann boundary (top, left) conditions
% u(x,y) = g(x,y) on a square grid using the Finite Difference method 5 point stencil.
%
% UC Berkeley Math 228B, Tyler Maltba (tyler_maltba@berkeley.edu)

% Vector of different grid sizes
nvec = [10 20 40];
ln = length(nvec);

% Allocation
Uerr = 0*(1:ln);
Qerr = Uerr;

% Parameters
L = 3.0;
H = 1.0;

% Vector of different choices of B
Bvec = [0, 0.5, 1.0];
lb = length(Bvec);

% Calculate errors of solution u and integral Q in reference domain 
for l = 1:lb
    B = Bvec(l);
    
    A = sqrt(0.25*(L-B)^2 - H^2);
    
    for k = 1:ln

        n = nvec(k);    
        [u_err, Q_err] = errPoisson(n,L,B,H); 
        Uerr(k) = u_err;
        Qerr(k) = Q_err;
        
        % plot solution u(x,y) for n = 40, B=0.5
        if B==0.5 && n==40
            
            [Q, u, xi, eta] = solvePoisson(n,L,B,H);
            
            % map back to physical domain
            x = (A*eta + B/2).*xi;
            y = H*eta;
            
            % plot solution u profile
            figure(2)
            subplot(1,2,1)
            contourf(x, y, u, 10);
            axis("equal"); colorbar();
            title("Contour: u(x,y)");
            set(gca, 'fontsize', 20, 'linewidth',1.75)
            
            subplot(1,2,2)
            surf(x, y, u);
            axis("equal"); colorbar();
            title("u(x,y)");
            shading interp
            set(gca, 'fontsize', 20, 'linewidth',1.75)
            shg
        end
    end
    
    figure(1)
    subplot(1,lb,l)
    loglog(1./nvec, 1./nvec.^2, '-b', 1./nvec, Uerr,'-.r',...
        1./nvec,Qerr,'--k','linewidth',2)
    title(['B = ', num2str(B)])
    set(gca, 'fontsize', 20, 'linewidth',1.75)
    xlabel('h')
    if l==1
        legend('1/h^2','E_u(h)','E_Q(h)','location','se')
    end
    shg
end

function [Q, u, xi, eta] = solvePoisson(n,L,B,H) 
% function for assemble the system matrix and RHS.
% Finite difference 9 point stencil.
% u is solution in reference domain, and Q is associated integral
    h = 1.0 / n;
    N = (n+1)^2;
    xi = h * (0:n);
    eta = xi;
    
    A = sqrt(0.25*(L-B)^2 - H^2);
    
    umap = reshape(1:N, n+1, n+1)';
    M = zeros(N, N);
    b = zeros(N, 1);
    
    for j = 1:n+1
        for i = 1:n+1
            
          % Solve -(D1*u_xi,xi + D2*u_eta,eta + D3*u_xi,eta + D4*u_xi) = 1
            
            D1 = (H^2 + (A*xi(i))^2) / (H^2*(A*eta(j) + B/2)^2);
            D2 = 1/H^2;
            D3 = -2*A*xi(i) / (H^2*(A*eta(j) + B/2));
            D4 = 2*xi(i)*A^2 / (H^2*(A*eta(j) + B/2)^2);
            
            row = umap(i,j);
            
            % left BC (Neumann, 2nd order, one-sided), u_xi = 0
            % -1/h * (1.5*u(1,j) - 2*u(2,j) + 0.5*u(3,j)) = 0 
            if i==1 && j<=n+1 && j>=2   
                M(row, row) = 1.5;
                M(row, umap(i+1,j)) = -2.0;     
                M(row, umap(i+2,j)) = 0.5;
                b(row) = 0;
                
            % right BC (Dirichlet), u = 0    
            elseif i==n+1 && j<=n+1 && j>=2     % right BC
                M(row, row) = 1.0;
                b(row) = 0;
            
             % bottom BC (Dirichlet), u = 0
            elseif j==1 && i<=n+1 && i>=1   
                M(row, row) = 1.0;
                b(row) = 0;    
             
            % top BC (Mixed Neumann, 2nd order), 
            % (A*eta + B/2)*u_eta - A*xi*u_xi = 0
            % central difference in xi direction
            % one-sided in eta direction
            % (A*eta + B/2)*(1.5*u(i,n+1) - 2*u(i,n) + 0.5*u(i,n-1))/h ...
            % -A*xi*(u(i+1,n+1)-u(i-1,n+1))/2h = 0 
            elseif j==n+1 && i<=n && i>=2  
                M(row, row) = 1.5*(A*eta(j) + B/2);
                M(row, umap(i,j-1)) = -2.0*(A*eta(j) + B/2);
                M(row, umap(i,j-2)) = 0.5*(A*eta(j) + B/2);
                M(row, umap(i+1,j)) = -A*xi(i)/2;
                M(row, umap(i-1,j)) = A*xi(i)/2;
                b(row) = 0;
                
            % 9-point stencil, 2nd order central differencing     
            else
                M(row, row) = 2*(D1+D2)/(h^2);
                M(row, umap(i+1,j)) = -(D1/(h^2) + D4/(2*h));
                M(row, umap(i-1,j)) = -(D1/(h^2) - D4/(2*h));
                M(row, umap(i,j+1)) = -D2/(h^2);
                M(row, umap(i,j-1)) = -D2/(h^2);                
                M(row, umap(i+1,j+1)) = -D3/(4*h^2);
                M(row, umap(i+1,j-1)) = D3/(4*h^2);
                M(row, umap(i-1,j+1)) = D3/(4*h^2);
                M(row, umap(i-1,j-1)) = -D3/(4*h^2);
                b(row) = 1.0;
            end
        end
    end
                
   %Sparse storage in Matlab
   M=sparse(M); 
   %Notice that matrix A depends the ordering of the global index, 
   %therefore it is not neccessary symmetric. 
   
   % numerical solution u on unit square
   u = reshape(M\b, n+1, n+1);
   
   % Find flow rate in reference domain   
   % Change of variables: |J(phi)|(xi,eta) = H*(A*eta + B/2)
   % int_phi(R) u(x,y) dxdy = int_R H*(A*eta + B/2)*u(xi,eta) d_xi d_eta
   % trapezoidal intergration
   
   eta1 = meshgrid(eta);
   Q = trapz(eta,trapz(xi,u*H.*(A*eta1 + B/2),2));
   
   % meshgrid of reference domain
   [xi,eta] = meshgrid(xi,eta);

end

function [u_err, Q_err] = errPoisson(n,L,B,H)
  
    % Numerical solution
    [Q, u, xi, eta] = solvePoisson(n,L,B,H); 
    
    % "true" solution with n = 80
    n0=80;
    [Q0, u0, xi0, eta0] = solvePoisson(n0,L,B,H);     
    
    % find u0 on the grid for u
    u0 = interp2(xi0,eta0,u0,xi,eta,'makima');
    
    % compute infinity norm of Error for u 
    u_err = max(max(abs(u-u0)));
    
    % Error for integral Q
    Q_err = abs(Q - Q0);
    
end
    
    
    
    
    
    
    
    