% Solve 2D Poisson's equation -(u_xx+u_yy) = f, with Dirichlet boundary condition
% u(x,y) = g(x,y) on a square grid using the Finite Difference method 5 point stencil.
%
% UC Berkeley Math 228B, Suncica Canic (canics@berkeley.edu)

[u, error] = testPoisson(100); % Change the resolution of the mesh
disp("Error infinity norm is: %d " + error);

function [A, b, x, y] = assemblePoisson(n, f, g) 
% function for assemble the system matrix and RHS.
% Finite difference 5 point stencil.
    h = 1.0 / n;
    N = (n+1)^2;
    x = h * (0:n);
    y = x;
    
    umap = reshape(1:N, n+1, n+1);
    A = zeros(N, N);
    b = zeros(N, 1);
    
    for j = 1:n+1
        for i = 1:n+1
            row = umap(i,j);
            if i==1 || i==n+1 || j==1 || j==n+1
                A(row, row) = 1.0;
                b(row) = g(x(i), y(j));
            else
                A(row, row) = 4.0;
                A(row, umap(i+1,j)) = -1.0;
                A(row, umap(i-1,j)) = -1.0;
                A(row, umap(i,j+1)) = -1.0;
                A(row, umap(i,j-1)) = -1.0;
                b(row) = f(x(i), y(j)) * h^2;
            end
        end
    end
                
   %Sparse storage in Matlab
   A=sparse(A); 
   %Notice that matrix A depends the ordering of the global index, 
   %therefore it is not neccessary symmetric. 
   

end
    
function [u, error] = testPoisson(n)

    uexact = @(x, y)exp(-4*(x-0.3).^2 - 9*(y-0.6).^2);
    f = @(x, y)exp(-4*(x-0.3).^2 - 9*(y-0.6).^2) * (26 - (18*y - 10.8).^2 - (8*x - 2.4).^2);
    [A, b, x, y] = assemblePoisson(n, f, uexact);
    
    u = reshape(A\b, n+1, n+1); % Numerical solution
    u0 = uexact(x',y);% Exact solution
    
    % subrountine for plotting
    subplot(2,2,1);
    contour(x, y, u, 10);
    contourf(x, y, u, 10);
    axis("equal"); colorbar();
    title("Contour: u_{computed}");
    
    subplot(2,2,2);
    contour(x, y, u0, 10);
    contourf(x, y, u0, 10);
    axis("equal"); colorbar();
    title("Contour: u_{exact}");
    
    subplot(2,2,3);
    surf(x, y, u);
    surf(x, y, u);
    axis("equal"); colorbar();
    title("u_{computed}");
    
    subplot(2,2,4);
    surf(x, y, u0);
    surf(x, y, u0);
    axis("equal"); colorbar();
    title("u_{exact}");
    
    %compute infinite norm
    error = max(max(abs(u-u0)));

end


    
    
    