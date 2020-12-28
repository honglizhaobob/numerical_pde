% Implementation of finite element method solving the 
% 2D Poisson problem on general deformed meshes

% with homogeneous neumann condition, Dirichlet BCs and 
% mixed conditions

% Hongli Zhao, honglizhaobob@berkeley.edu
function u = fempoi(p, t, e)

% Main routine to compute the numerical solutions for the Poisson
% 2D problem, using finite element analysis. The script first
% assembles a global stiffness matrix using general neumann condition
%                       n * nabla u = 0
% and replaces the entries for Dirichlet if e is nonempty
%
%   Inputs:
%
%       p = nx2, coordinates of nodes
%       t = Nx3, T matrix (associates each triangle with 3 
%           numbered vertices)
%       e = an array of nodes that have homogenous Dirichlet BC

%   Outputs:
%
%       u = nxn matrix with data representing the numerical surface

% number of nodes
M = size(p,1);

% number of triangles
% note: M < N since one node can be used 
% possibly for many triangles
N = size(t,1);


% size of Dirichlet boundary
G = size(e);

A = zeros(M,M); b = zeros(M,1);
for n = 1:N
    [An, bn] = elem_stiff(p,t,n);
    for alpha = 1:3
        for beta = 1:3
            % set entries
            if ~ismember(t(n,alpha),e)                
                A(t(n,alpha), t(n,beta)) =...
                    A(t(n,alpha), t(n,beta)) + An(alpha,beta);       
            end
        end
        b(t(n,:)) = b(t(n,:)) + bn(alpha);
    end
end

for i = 1:G
    replace=e(i);
    % any nonzero would work
    A(replace,replace) = 1;
end

A = sparse(A);

xi = A\b;
u = xi;

end

% ==============================
%          Subroutines
% ==============================
function tri = query(p,t,num)
        % returns the 3 (xi,yi) coords of a triangle
        % 3x2 matrix
        tri = p(t(num, :), :);
end

function [An, bn] = elem_stiff(p,t,n)
    tri = query(p,t,n);
    B = [tri ones(3,1)];
    % area of triangle
    area_K = (1/2)*det(B);
    C = inv(B);
    An = zeros(3);
    for i = 1:3
        for j = 1:3
            grad = sum(C(1:2,i).*C(1:2,j)); 
            An(i,j) = grad * area_K;
        end
    end 
    bn = ones(3,1) * (1/6)*det(B);
end