% contains code to simulate the Isentropic Vortex Convection problem,
% using our Pade solver and compact filter.
%
clear
% control sequences
with_plot = 0;
% final time

% change N for different results
N = 128; % N = 32,64,128
h = 10.0/N;
x = h * (0:N);
y = x;
[x,y] = meshgrid(x,y);
T = 5 * sqrt(2);
k = 0.3*h;
adjust_factor = ceil(T/k); % factor to adjust k
k = T/adjust_factor;
%===================== parameters
x_c = 5; y_c = 5;
r_inf = 1; u_inf = 0.1; 
v_inf = 0; p_inf = 1;
gamma = 7/5; b = 0.5;

alpha = 0.499; % 0.48, 0.499

FinalTime = T;

rr = sqrt((x - x_c).^2 + (y - y_c).^2);
r = (1 - ((b^2)*(gamma-1)/(8*gamma*pi^2))*exp(1-rr.^2)).^(1/(gamma-1));
u = 0.1 - (0.5*b/pi)*exp(0.5*(1-rr.^2)).*(y-y_c);
v = (0.5*b/pi)*exp(0.5*(1-rr.^2)).*(x-x_c);
p = r.^gamma;

ru = r.*u;
rv = r.*v;
rE = p/(gamma-1) + r.*(u.^2 + v.^2)/2;
%===================== function calls



% with all other dependencies, integrate over time to T
for t = 1:adjust_factor
    [r,ru,rv,rE] = euler_rk4step(r,ru,rv,rE,h,k,alpha);
end

if with_plot ~= 0
    % see what happened
    figure(1)
    subplot(2,2,1)
    contourf(x,y,r,16)
    shading interp
    title('\rho')
    colorbar

    subplot(2,2,2)
    contourf(x,y,ru,16)
    shading interp
    title('\rho u')
    colorbar

    subplot(2,2,3)
    contourf(x,y,rv,16)
    shading interp
    title('\rho v')
    colorbar

    subplot(2,2,4)
    contourf(x,y,rE,16)
    shading interp
    title('\rho E')
    colorbar
end


% errors visualization
x_exact = x_c + T*u_inf;
d = sqrt((x-x_exact).^2 + (y-y_c).^2);
r_exact = (1 - ((b^2)*(gamma-1)/(8*gamma*pi^2))*exp(1-d.^2)).^(1/(gamma-1));
u_exact = u_inf - (0.5*b/pi)*exp(0.5*(1-d.^2)).*(y-y_c);
v_exact = v_inf + (0.5*b/pi)*exp(0.5*(1-d.^2)).*(x-x_exact);
p_exact = r_exact.^gamma;

ru_exact = r_exact.*u_exact;
rv_exact = r_exact.*v_exact;
rE_exact = p_exact/(gamma-1) + r_exact.*(u_exact.^2 + v_exact.^2)/2;
% generate exact solutions

figure(2)
subplot(2,2,1)
contourf(x,y,r_exact,16)
shading interp
title('\rho')
colorbar

subplot(2,2,2)
contourf(x,y,ru_exact,16)
shading interp
title('\rho u')
colorbar

subplot(2,2,3)
contourf(x,y,rv_exact,16)
shading interp
title('\rho v')
colorbar

subplot(2,2,4)
contourf(x,y,rE_exact,16)
shading interp
title('\rho E')
colorbar


r_error_mat = r - r_exact;
ru_error_mat = ru - ru_exact;
rv_error_mat = rv - rv_exact;
rE_error_mat = rE - rE_exact;

% calculate infinity norms for errors
r_error = max(max(abs(r_error_mat)));
ru_error = max(max(abs(ru_error_mat)))
rv_error = max(max(abs(rv_error_mat)))
rE_error = max(max(abs(rE_error_mat)))

disp(['rho error is: ', num2str(r_error)]);
disp(['rho u error is: ', num2str(ru_error)]);
disp(['rho v error is: ', num2str(rv_error)]);
disp(['rho E error is: ', num2str(rE_error)]);