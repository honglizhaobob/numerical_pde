% Driver script for solving the 1D Burgers equations using linear basis
% functions.
% Both strong form and weak form are implemented.
% Option=1: strong form, 2: weakform. 
clc; clear; clf;
global N J nx K M lmln 
global vmapM vmapP vmapI vmapO


N = 10; %Number of elements 4,10,20,50
xL=-1; xR=1;
X = linspace(xL, xR, N+1);
FinalTime =0.5;

% We use linear basis function.
% We compute the mass matrix on the standard/reference element [0,1]; 
% All elements are mapped to this standard element;
syms t;
phi=[(1-t)/2, (1+t)/2];
phip=[diff(phi(1),1);diff(phi(2),1)];

J =((xR-xL)/N/2);
M=zeros(2,2);
for i=1:2
    for j=1:2
   M(i,j)=int(phi(i)*phi(j)*J,t,-1,1);     
    end
end

K=zeros(2,2);
for i=1:2
    for j=1:2
   K(i,j)=int(phi(i)*phip(j),t,-1,1);     
    end
end


r=[-1;1];
x = ones(2,1)*X(1:N) + 0.5*(r+1)*(X(2:N+1)-X(1:N));%find the x coordinates for all the elements


disp('=> x vector is: ')
disp((x))

disp('=> mass matrix M')
disp(M)

disp('=> matrix K')
disp(K)




%Flux lifting term, i.e., the boundary term in weak formulation
lmln = zeros(2,2); 
lmln(1,1) = 1.0; lmln(2,2) = 1.0;

% Surface normals for each element, direction points to the right are positive
nx = zeros(2, N); nx(1, :) = -1.0; nx(2, :) = 1.0;

%Connectivity
vmapM   = zeros(2*N,1);  
vmapP   = zeros(2*N,1); 
for i=1:1:N
vmapM(i*2-1,1)=(i-1)*2+1;
vmapM(i*2,1)=(i-1)*2+2;
end
for i=1:1:N
vmapP(i*2-1,1)=(i-1)*2;
vmapP(i*2,1)=(i-1)*2+2+1;
end
vmapP(1,1)=1; vmapP(2*N,1)=N*2;
vmapI = 1; vmapO = 2*N;% Create specific left (inflow) and right (outflow) maps





% Initial condition plotting
u = sin(pi*x);

figure(1);
plot(x,u,'r'); 
hold on;

disp('=> initial u')
disp(u)


time = 0;
% compute time step size
xmin = min(abs(x(1,:)-x(2,:)));
CFL=0.2; umax = max(max(abs(u)));
dt = CFL* min(xmin/umax);
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps; 
% for plotting
disp('=> time step dt is: ')
disp(dt)

disp('=> *** size of x')
disp(size(x))
disp('=> *** size of u')
disp(size(u))
plot_vector = [];

for tstep=1:Nsteps
    rhs = BurgersRHS(u); 
    
    if time == 0
        figure(2); hold on; grid on;
        plot(x,u, 'b','LineWidth',2);
        axis([-1 1 -1.5 1.5]);
        title(['N = ', num2str(N), ' solution u at various times'])
        xlabel('x axis');
        ylabel('u axis');
        disp('=> plotted')
        disp(['=>=> plotted time is: ', num2str(time)])
    end
    if time == 0.25 || tstep == ceil(Nsteps/2)
       figure(2); hold on; grid on;
       plot(x, u, 'm','LineWidth',2);
    end
    if tstep == Nsteps
       figure(2); hold on; grid on
       plot(x,u, 'black','LineWidth',2);
    end
    
    
    
    
    
    
    
    
    
    u = u+dt*rhs; %Forward Euler to update u^{n+1} from u^n
    time = time+dt;
    %if (mod(tstep,10))==1 %Plotting every 10th time step
    %if time == 0.25 || tstep == Nsteps % plotting at only 3 times
       %plot(x,u,'b-');
       %axis([-1 1 -1.5 1.5]);
       %title(['N = ', num2str(N), ' Solution at time t = ',num2str(time)]);
       %title(['N = ', num2str(N), ' solution u at various times'])
       %xlabel('x axis');
       %ylabel('u axis');
       %disp('=> plotted')
       %disp(['=>=> plotted time is: ', num2str(time)])
    %pause
    %end
end


function [rhs] = BurgersRHS(u)
global N nx K M lmln 
global vmapM vmapP vmapI vmapO

 du = zeros(2,N); du(:) = u(vmapM)-u(vmapP);
 MaxVel = max(max(abs(u)));% Wave speed=constant C in Lax-Friedrichs flux
 
     %Right hand-side of the weak form
     % TODO: IMPLEMENT THE LAX-FRIEDRICHS FLUX
     rhs = zeros(2,N);
     flux = rhs;
     
     % fixme
     c = MaxVel;
     % u0 and u1's at each element
     u0 = u(1,:); u1 = u(2,:);
     for i = 1:N % loop over each element
         F = zeros(2,1);
         if i == 1 % at left boundary
             u1km1 = u0(1); 
             u0k = u0(1);           
             f1 = (1/4)*(u1km1^2 + u0k^2) + (c/2)*(u1km1-u0k);
             u1k = u1(1);
             u0kp1 = u0(2);
             f2 = (1/4)*(u1k^2 + u0kp1^2) + (c/2)*(u1k-u0kp1);
             F(1) = f1; F(2) = -f2;
         elseif i == N % at right boundary
             u1km1 = u1(end-1);
             u0k = u0(end);
             f1 = (1/4)*(u1km1^2 + u0k^2) + (c/2)*(u1km1-u0k);
             u1k = u1(end);
             u0kp1 = u1k;
             f2 = (1/4)*(u1k^2 + u0kp1^2) + (c/2)*(u1k-u0kp1);
             F(1) = f1; F(2) = -f2;
         else
             % else points in the middle
             u1km1 = u1(i-1);
             u0k = u0(i);
             f1 = (1/4)*(u1km1^2 + u0k^2) + (c/2)*(u1km1-u0k);
             u1k = u1(i);
             u0kp1 = u0(i+1);
             f2 = (1/4)*(u1k^2 + u0kp1^2) + (c/2)*(u1k-u0kp1);
             F(1) = f1; F(2) = -f2;
         end
         flux(:,i) = F;
     end
     % with the lax friedrichs part, assemble 
     % rhs
     for i = 1:N % loop over each elem
         rhs(:,i) = 0.5*inv(M)*transpose(K)*(u(:,i).^2) + ...
             inv(M)*flux(:,i);
     end
end




