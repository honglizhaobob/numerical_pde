% u''''(x) = f(x) problem using FEM
% compare exact solution with numerical solution

% dictionary
syms phi1k1(x) phi1k2(x) phi2k1(x) phi2k2(x) f(x)
% right hand side
f(x) = 480.*x - 120;
% phi1 first half + second half
phi1k1(x) = 4.*x.^3-2.*x.^2;
phi1k2(x) = 4.*x.^3 - 10 .* x.^2 + 8.*x - 2;
% phi2 first half + second half
phi2k1(x) = -12*x.^3 + 10.*x.^2;
phi2k2(x) = 20.*x.^3 - 46.*x.^2 + 32.*x - 6;

% dictionary for second ord derivatives
p11pp = diff(diff(phi1k1(x)));
p12pp = diff(diff(phi1k2(x)));
p21pp = diff(diff(phi2k1(x)));
p22pp = diff(diff(phi2k2(x)));

% stiff a11, convert to function handle first
% phi1'' * phi1''
a11_firsthalf = matlabFunction(p11pp * p11pp);
a11_secondhalf = matlabFunction(p12pp * p12pp);
a11 = integral(a11_firsthalf, 0, 1/2) + ...
    integral(a11_secondhalf, 1/2, 1);

% stiff a12
% phi1'' * phi2''
a12_firsthalf = matlabFunction(p11pp * p21pp);
a12_secondhalf = matlabFunction(p12pp * p22pp);
a12 = integral(a12_firsthalf, 0, 1/2) + ...
    integral(a12_secondhalf, 1/2, 1);

% stiff a21
% phi2'' * phi1''

a21 = a12; % symmetric

% stiff a22
% phi2'' * phi2''
a22_firsthalf = matlabFunction(p21pp * p21pp);
a22_secondhalf = matlabFunction(p22pp * p22pp);
a22 = integral(a22_firsthalf, 0, 1/2) + ...
    integral(a22_secondhalf, 1/2, 1);

% stiffness matrix
A = [a11,a12;a21,a22];

% now load vector
b1_firsthalf = matlabFunction(phi1k1(x) * f(x));
b1_secondhalf = matlabFunction(phi1k2(x) * f(x));
b1 = integral(b1_firsthalf, 0, 1/2) + ...
    integral(b1_secondhalf, 1/2, 1);

% b2
b2_firsthalf = matlabFunction(phi2k1(x) * f(x));
b2_secondhalf = matlabFunction(phi2k2(x) * f(x));
b2 = integral(b2_firsthalf, 0, 1/2) + ...
    integral(b2_secondhalf, 1/2, 1);

% load vector b
b = [b1;b2];

% find xi
xi = A\b;

% with xi, can get function
% phi1 and phi2 are piecewise P3 functions
phi1 = piecewise(0<=x<=0.5, phi1k1, 0.5<x<=1, phi1k2);
phi2 = piecewise(0<=x<=0.5, phi2k1, 0.5<x<=1, phi2k2);

% uh = xi(1)*phi1 + xi(2)*phi2
uh = xi(1) * phi1 + xi(2) * phi2;

% now begin plotting solutions
syms uexact(x)
uexact(x) = 4.*x.^5-5.*x.^4-2.*x.^3+3.*x.^2;

% generate mesh for visualization
xx = linspace(0,1,100);
uh_data = double(uh(xx));
uexact_data = double(uexact(xx));

figure(1)

plot(xx,uh_data,'r',xx,uexact_data,'b','LineWidth',1.2)
grid on
title('FEM Numerical Solution vs Exact Solution')
xlabel('x'); ylabel('u'); 
legend('u_h', 'u_{exact}')