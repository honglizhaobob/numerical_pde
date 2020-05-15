% routine to calculate finest level solution
% and compute error vectors
function errors=poiconv(pv,hmax,nrefmax)
    % use the finest level as exact solution
    [p_exact,t_exact,e_exact] = pmesh(pv, hmax, nrefmax);
    % solve on finest level using fempoi first
    U_exact = fempoi(p_exact,t_exact,e_exact);
    % allocate error vector and iteratively update
    errors = zeros(nrefmax,1);
    for i = 0:nrefmax-1
        [p_test,t_test,e_test] = pmesh(pv,hmax,i);
        U_test = fempoi(p_test,t_test,e_test);
        % take the first n solutions 
        errors(i+1) = max(abs(U_test(1:size(p_test,1))-...
            U_exact(1:size(p_test,1))));
    end  
end