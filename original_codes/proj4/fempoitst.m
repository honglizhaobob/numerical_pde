% Square, Dirichlet left/bottom
pv = [0,0; 1,0; 1,1; 0,1; 0,0];
[p,t,e] = pmesh(pv, 0.1, 0);
e = e(p(e,1) < 1e-6 | p(e,2) < 1e-6);
u = fempoi(p,t,e);

subplot(2,2,1);
trisurf(t,p(:,1),p(:,2),u);

% Circle, all Dirichlet
n = 32; phi = 2*pi*(0:n)'/n;
pv = [cos(phi), sin(phi)];
[p,t,e] = pmesh(pv, 2*pi/n, 0);
u = fempoi(p,t,e);

subplot(2,2,2);
trisurf(t,p(:,1),p(:,2),u);

% Generic polygon geometry, mixed Dirichlet/Neumann
x = (0:.1:1)';
y = 0.1 * cos(10*pi*x);
pv = [x,y; .5,.6; 0,.1];
[p,t,e] = pmesh(pv, 0.04, 0);
e = e(p(e,2) > .6 - abs(p(e,1) - 0.5) - 1e-6);
u = fempoi(p,t,e);

subplot(2,2,3);
trisurf(t,p(:,1),p(:,2),u);

