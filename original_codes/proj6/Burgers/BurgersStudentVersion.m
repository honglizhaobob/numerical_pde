% Driver script for solving the 1D Burgers equations
%clc; clear; clf;
global N P
global V inV Dr LIFT rx J nx Fscale
global vmapM vmapP vmapI vmapO

% 1-20, 3-20, 1-50, 3-50

figure(1); hold on; grid on;
pp = [1,3];
nn = [20,50];
for ii = 1:2
    for jj = 1:2
%=== added for loop

P = pp(ii);%Polynomial order
N = nn(jj); %Number of elements
if P == 3 && N == 50
    continue;
end

% added for loop

xL=-1; xR=1;
FinalTime =0.5;

X = linspace(xL, xR, N+1);
Np = P+1;

% Compute basic Legendre Gauss Lobatto grid
r = JacobiGL(0,0,P);

% Build reference element matrices
[V,inV,Dr]  = VandermondeMatrix(P, r); 

lmln = zeros(Np,2);%l_m()*l_n()
lmln(1,1) = 1.0; lmln(Np,2) = 1.0;

% Lift operator: surface integral 
% corresponding to Minverse times surface terms
LIFT = V*(V'*lmln);

%find the x coordinates for all the elements
x = ones(Np,1)*X(1:N) + 0.5*(r+1)*(X(2:N+1)-X(1:N));

Edge_Idx  = [1, Np];% Elemental edge node index

% Surface normals for each element, 
% direction points to the right is positive
nx = zeros(2, N); nx(1, :) = -1.0; nx(2, :) = 1.0;
J = Dr*x; rx = 1./J; % Jacobian
Fscale = 1./(J(Edge_Idx,:));

%Connectivity
vmapM   = zeros(2*N,1); 
vmapP   = zeros(2*N,1); 
for i=1:1:N
vmapM(i*2-1,1)=(i-1)*Np+1;
vmapM(i*2,1)=(i-1)*Np+Np;
end
for i=1:1:N
vmapP(i*2-1,1)=(i-1)*Np;
vmapP(i*2,1)=(i-1)*Np+Np+1;
end
vmapP(1,1)=1;
vmapP(2*N,1)=N*Np;

% Create specific left (inflow) and right (outflow) maps
vmapI = 1; vmapO = N*Np;

% Initial conditions sin(pi*x)
u = sin(pi*x);
figure(1);  %plot(x,u,'LineWidth',1.4);hold on;

time = 0;

% compute time step size
xmin = min(abs(x(1,:)-x(2,:)));
CFL=0.25; umax = max(max(abs(u)));
dt = CFL* min(xmin/umax);
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps; 


for tstep=1:Nsteps
    rhs = BurgersRHS(u); 
    u = u+dt*rhs; %Forward Euler to update u^{n+1} from u^n
    u  = Limiter(u);
    time = time+dt;
    %if (mod(tstep,10))==1
    if tstep == Nsteps % only plot at certain times
        disp('=> plotting')
        disp(['current time: ', num2str(time)])
        if tstep == ceil(Nsteps/2)
            if P == 3 && N == 20
                disp(['==> plotting: '...
                    , num2str(P), '++', num2str(N)])
                plot(x(:),u(:),'b','LineWidth',1.5); hold on; grid on;
            elseif P == 1 && N == 20
                disp(['==> plotting: '...
                    , num2str(P), '++', num2str(N)])
                plot(x(:),u(:),'m','LineWidth',1.5); hold on; grid on;
            else
                disp(['==> plotting: '...
                    , num2str(P), '++', num2str(N)]);
                plot(x(:),u(:),'black','LineWidth',1.5); hold on; grid on;
            end
        else
            %plot(x(:),u(:),'black','LineWidth',1.5); hold on; grid on;
            if P == 3 && N == 20
                disp(['==> plotting: '...
                    , num2str(P), '++', num2str(N)])
                plot(x(:),u(:),'b','LineWidth',1.5); hold on; grid on;
            elseif P == 1 && N == 20
                disp(['==> plotting: '...
                    , num2str(P), '++', num2str(N)])
                plot(x(:),u(:),'m','LineWidth',1.5); hold on; grid on;
            else
                disp(['==> plotting: '...
                    , num2str(P), '++', num2str(N)]);
                plot(x(:),u(:),'black','LineWidth',1.5); hold on; grid on;
            end
        end
        title('u(x,0.5) for diff methods');
        xlabel('x axis');
        ylabel('u axis');
    end
end

%=== added for loop
    end
    
end
legend('P = 1, N = 20','P = 1, N = 50','P = 3, N = 20','Location','nw')




function [rhs] = BurgersRHS(u)
global N Dr LIFT rx nx Fscale
global vmapM vmapP vmapI vmapO

% Compute Jump in u
du = zeros(2,N); 
du(:) = u(vmapM)-u(vmapP);

% Evaluate flux (u^2)/2
du2 = zeros(2,N); du2(:) = (u(vmapM).^2-u(vmapP).^2)/2.0;

MaxVel = max(max(abs(u)));% Wave speed

% Lax-Friedrichs numerical flux 
flux = nx.*(du2/2.0) - MaxVel/2.0.*du ;

% Compute RHS
rhs = -(rx.*(Dr*(u.^2/2)) - LIFT*(Fscale.*flux));
return
end

function ulimit = Limiter(u) %Slope Limiter
% Purpose: Apply slopelimiter (Pi^N) 
% to u assuming u an N'th order polynomial         
global N P V inV 

% Compute cell averages
uh = inV*u; uh(2:P+1,:)=0; uavg = V*uh; v = uavg(1,:);

% Apply slope limiter as needed.
ulimit = u; eps0=1.0e-8;

% find end values of each element
ue1 = u(1,:); ue2 = u(end,:);

% find cell averages
vk = v; vkm1 = [v(1),v(1:N-1)]; vkp1 = [v(2:N),v(N)]; 

% Apply reconstruction to find elements in need of limiting
ve1 = vk - minmod([(vk-ue1);vk-vkm1;vkp1-vk]);
ve2 = vk + minmod([(ue2-vk);vk-vkm1;vkp1-vk]);
ids = find(abs(ve1-ue1)>eps0 | abs(ve2-ue2)>eps0);
end
