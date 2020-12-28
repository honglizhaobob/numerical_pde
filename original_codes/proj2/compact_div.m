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