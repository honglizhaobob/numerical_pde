function u = fempoi2(p,t,e)
    ntri = size(t,1);
    np = size(p,1);
    A = [];
    b = zeros(np,1);
    tmp = 0;

    for i = 1:ntri
        V=[ones(3,1), p(t(i, :), :)];
        area = (1/2)*det(V);
        C = inv(V);
        for j = 1:3
            for k = 1:3
                tmp = area*(sum(C(2:3,j).*C(2:3,k)));
                if ~ismember(t(i,j),e)
                    A=[A;t(i, j), t(i, k), tmp];
                end
            end
            b(t(i)) = b(t(i)) + (1/6)*det(V);
        end
        
    end
    tt=size(e);
    for i = 1:tt
        %disp(num2str(tmp))
        A=[A;e(i),e(i),1];
        b(i) = 0.0;
    end
    A
    s = size(A,1);
    matrix = zeros(np,np);
    for i=1:s
       x=A(i,1);
       y=A(i,2);
       matrix(x,y)=matrix(x,y)+A(i,3);
    end
    u=matrix\b;
end
