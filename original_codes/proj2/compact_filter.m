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
end