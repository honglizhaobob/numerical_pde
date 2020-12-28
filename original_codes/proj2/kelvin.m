% problem 2, kelvin-helmholtz instability on 
% unit domain
N = 256;
h = 1.0/N;

x = (0:h:1.0); y = (0:h:1.0);
T = 1.0;
k = 0.3*h;
adjust_factor = ceil(T/k); % factor to adjust k
k = T/adjust_factor;

[x,y] = meshgrid(x,y);
g = 7/5; % gas constant
alpha = 0.48; % 0.48, 0.499
% initial conditions
r = zeros(size(x));
F1 = find(abs(y-0.5) < (0.15 + sin(2*pi*x)/200));
F2 = find(~(abs(y-0.5)<0.15 + sin(2*pi*x)/200));
r(F1) = 2;
r(F2) = 1;
u = r - 1;
v = zeros(size(r));
p = 3*ones(size(r));

ru = r.*u;
rv = r.*v;
rE = p/(g-1) + r.*(u.^2 + v.^2)/2;

% integrate in time
for t = 0:k:T
    [r,ru,rv,rE] = euler_rk4step(r,ru,rv,rE,h,k,alpha); 
end
    figure(2)
    contourf(x,y,r,16)
    shading interp
    title('\rho')
    colorbar

