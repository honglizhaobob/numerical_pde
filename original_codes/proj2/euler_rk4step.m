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