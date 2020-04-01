%%%%%%%%%%%%%%%%%%%%%%%%%
%        2 steps Lax-Wendroff scheme for              %
%                    1D AQ system                           %
%                        4/3/2019                              %
%%%%%%%%%%%%%%%%%%%%%%%%%

global dx alpha nu rho sigma E R0 h G0 A0 L c dt T1 Xin Xout

figure(1); clf;
figure(2); clf;

% Read input data obtained using fft
Xin = dlmread('Hyperaemia_ProximalPressure_fft.dat');
Xout = dlmread('Hyperaemia_DistalPressure_fft.dat');
externalPressure =(75.21+(179.48-75.21)/3)*1333;%(Max-Min)/3 in dynes/cm^2 
disp(externalPressure)

N=26;%# of elements
alpha = 4/3;
nu = 0.035; %viscosity (cm^2/s)
rho = 1.055; %density (g/cm^3)
sigma = 0.49; %Poi sson ratio
E = 1e7; %4e6; %Young's modulus (dynes/cm^2)
R0 = 0.2; %Dilatated Radius (cm) due to administration of adenosine
h = 0.0475; %Thickness (cm)
G0 = (E*h)/R0/(1-sigma^2);
A0 = R0^2;
L = 4.5; %Artery length (cm)
dx= L/N;
CFLfactor=0.6;
c = sqrt(G0/rho); %sound speed
dt = dx/c*CFLfactor; %CFL
% Here we used only c for the CFL condition, instead of the full 
%eigenvalues lambda_1 or lambda_2 because the other terms in the 
% eigenvalue calculation are much smaller that c and do not
%influence the CFL condition much

Ncycles=3; %# of cycles we are calculating; 
% it takes at least one cycle
% for the solution to "stabilize" (transient period)

T1=.763; %Time interval for one cycle
T = T1*Ncycles;
x = 0:dx:L;
t = 0:dt:T;
index=1; %movie frame index
index2=1;


A = zeros(1,length(x))+R0^2; %A at time level n
Q = zeros(1,length(x)); %Q at time level n
A1= zeros(1,length(x)); %A at time level n+1
Q1= zeros(1,length(x)); %Q at time level n+1
U = zeros(2,length(x));
U1 = zeros(2,length(x));
U(1,:) = A;
U(2,:) = Q;
maxflow = Q(1,N/2); % records Q_max at x=N/2, update every time loop
% ==============================
%       Main loop solver
% ==============================
for n = 1:length(t)-1
    % report maxflow over one cycle (2nd cycle)
    if t(n) >= 1 * T1 && t(n) <= 2 * T1 && maxflow_init == 0
        disp('==========*==========*==========*==========')
        disp('started tracking maxflow at midpoint')
        maxflow = Q(1,N/2);
        disp(['current maxflow= ', num2str(maxflow)])
        maxflow_init = 1;
        % === reporting max Q(1,N/2)
        
    end
    if t(n) >= 1 * T1 && t(n) <= 2 * T1
        if Q(1,N/2) > maxflow
           maxflow = Q(1,N/2);
           disp("... reporting max flow at midpoint")
           disp(['current max flow = ', num2str(Q(1,N/2))]);
        end
    end
    for i = 1:length(x)
        if i == 1 %Left boundary (x=0) 
            
            % CALCULATE A1 and Q1 using upwinding;
            
            %U1(1,i) = A1(i);
            %U1(2,i) = Q1(i);
            
            % ==========
            %Write the code for boundary data at 0
            %U1(1,i) = A1(i);
            %U1(2,i) = Q1(i);
            
            % ==================== FIXME
            % begin upwind
            u_star = U(:,i+1); % this should be u_2^n
            u_0 = U(:,i); % known u_1^n, right at boundary x=0
            disp(['value of A_star is: ', num2str(u_star(1))])
            disp(['value of Q_star is: ', num2str(u_star(2))])
            % get jacobian matrix B for this U*
            B = jacobi(u_star);
            
            %=== added
            %B0 = jacobi(u_0);
            
            %===
         
            % get left eigenvectors and eigenvalues
            % eigenvalues not necessarily in ascending
            % order, need reordering
            [~,D,W] = eig(B);
            %=== new added
            %[~,D0,W0] = eig(B0);
            %eigvals0 = diag(D0); % unsorted
            %[d0,ind0] = sort(eigvals0,'ascend'); 
            %Ds0 = D0(ind0,ind0);
            %Ws0 = W0(:,ind0);
            %LM0 = Ws0';
            %===
            
            % get reordering
            [d,ind] = sort(diag(D),'ascend'); 
            Ds = D(ind,ind);
            Ws = W(:,ind);
            % left eigenvector matrix LM
            % LM * B = D * LM
            LM = Ws';
            lambdas = d;
            lambda1 = lambdas(1);
            lambda2 = lambdas(2);
            
            %[lambda2,lambda1] = lambda_pair(u_star);
            
            disp(['value of lambda1 is: ',...
                num2str(lambda1)])
            disp(['value of lambda2 is: ',...
                num2str(lambda2)])
            
            
            % begin main steps for upwind
            lu_star = LM * u_star; % [w,v]_2^n
            

            %=== new switch
            lu_0 = LM * u_0; % [w,v]_1^n
            %lu_0 = LM0 * u_0; % [w,v]_1^n
            %===
            
            % get appropriate source terms, 
            % see writeup for details
            
            %=== new switch
            src_0 = LM * rhs(u_0); 
            %src_0 = LM0 * rhs(u_0); 
            %===
            src_star = LM * rhs(u_star);
            
            % now use upwind to get [w,v]_1^{n+1}
            % -=============
            % need FIX
            w0 = lu_0(1);
            w1 = lu_star(1);
            
            s1_1 = src_star(1);
            s1_0 = src_0(1);
            
            w_boundary = ...
                w0 - (dt/dx)*lambda1*(w1-w0) +...
                (dt/2)*(s1_1+s1_0);
            
            % calc A precisely
            A_boundary = A0 * (Inlet(t(n))/G0 + 1)^2;
            % we've now recovered u at boundary
            A1(i) = A_boundary; %Q1(i) = Q_boundary;
            
            % w_boundary is a linear combination of
            % A_boundary and Q_boundary
            % with A_boundary, can calculate 
            % Q_boundary
            l_vec = LM(1,:);
            l_11 = l_vec(1); l_12 = l_vec(2);
            % l_11 * A + l_12 * Q = w
            % therefore: Q = (w - l_11*A)/l_12
            
            Q_boundary = (w_boundary - l_11*A_boundary)/l_12;
            Q1(i) = Q_boundary;
            
            
           
            %   Old: 
            
            %lu_boundary = ...
            %    lu_L - (dt/dx)*lambda1*(lu_star-lu_L)+...
            %    (dt/2)*(src_0 + src_star);
            
            
            % ==============
            
            % can now recover Q at the boundary, 
            % not A, because A can be calced
            % more precisely and directly from
            % the pressure formula at boundary
            
            %u_boundary = LM\lu_boundary;
            %Q_boundary = u_boundary(2);
            
            U1(:,i) = [A_boundary; Q_boundary];
            
        elseif x(i) == L  %Right boundary (x=L)
             
            % CALCULATE A1 and Q1 using upwinding;
              
            %U1(1,i) = A1(i);
            %U1(2,i) = Q1(i);
            
            % ==============
            %Write the code for boundary data at L
              
            % U1(1,i) = A1(i);
            % U1(2,i) = Q1(i);
            % ==================== FIXME
            % code remains largely similar
            % to the left boundary
            % begin upwind
            u_star = U(:,i-1); % this should be u_{Nx-1}^n
            u_L = U(:,i); % known u_Nx^n, right at boundary x=L
            disp(['    value of A_star is: ', num2str(u_star(1))])
            disp(['    value of Q_star is: ', num2str(u_star(2))])
            % get jacobian matrix B for this U*
            B = jacobi(u_star);
            
            %==== new added
            %Bn = jacobi(u_L);
            %[~,Dn,Wn] = eig(Bn);
            %eigvalsn = diag(Dn); % unsorted
            %[dn,indn] = sort(eigvalsn,'ascend'); 
            %Dsn = Dn(indn,indn);
            %Wsn = Wn(:,indn);
            %LMn = Wsn';
            %====
            
            % get left eigenvectors and eigenvalues
            % eigenvalues not necessarily in descending
            % order, need reordering
            [~,D,W] = eig(B);
            
            % get reordering
            [d,ind] = sort(diag(D),'ascend'); 
            Ds = D(ind,ind);
            Ws = W(:,ind);
            % left eigenvector matrix LM
            % LM * B = D * LM
            LM = Ws';
              
            lambdas = d;
            lambda1 = lambdas(1);
            lambda2 = lambdas(2);
            
            %[lambda2,lambda1] = lambda_pair(u_star);
            
            disp(['value of lambda1 is: ', num2str(lambda1)])
            disp(['value of lambda2 is: ', num2str(lambda2)])
            
            % begin main steps for upwind
            lu_star = LM * u_star; % [w,v]_2^n
            
            %=== new switch
            lu_L = LM * u_L; % [w,v]_1^n
            %lu_L = LMn * u_L; % [w,v]_1^n
            %===
            
            % get appropriate source terms, 
            % see writeup for details
            
            %=== new switch
            src_L = LM * rhs(u_L); 
            %src_L = LMn * rhs(u_L); 
            %===
            src_star = LM * rhs(u_star);
            
            % -=============
            % need fix
            vNm1 = lu_star(2);
            vN = lu_L(2);
            s2_Nm1 = src_star(2);
            s2_N = src_L(2);
            
            v_boundary = ...
                vN - (dt/dx)*lambda2*(vN-vNm1) +...
                (dt/2)*(s2_N + s2_Nm1);
            
            % calc A precisely
            A_boundary = A0 * (Outlet(t(n))/G0 + 1)^2;
            A1(i) = A_boundary; 
            
            l_vec = LM(2,:);
            % l_21 * A + l_22 * Q = v
            % Q = (v - l_21 * A)/l_22
            
            l_21 = l_vec(1); l_22 = l_vec(2);
            Q_boundary = (v_boundary - l_21*A_boundary)/l_22;
            Q1(i) = Q_boundary;
                      
            % we've now recovered u at boundary
            
            
            %   Old:
            % now use upwind to get [w,v]_1^{n+1}
            %lu_boundary = ...
            %    lu_star - (dt/dx)*lambda2*(lu_L - lu_star)+...
            %    (dt/2)*(src_L + src_star);
            
            % -============= old below
            
            % can now recover Q at the boundary, 
            % not A, because A can be calced
            % more precisely and directly from
            % the pressure formula at boundary
            
            %u_boundary = LM\lu_boundary;
            %Q_boundary = u_boundary(2);
            %Q1(i) = Q_boundary;
            
            % ====== old above
            U1(:,i) = [A_boundary; Q_boundary];
        
        else %Inner elements (using the two step lax-wendroff scheme with Strang splitting)
             
            %U_midA = U(:,i);  % MODIFY
            %U_midB = U(:,i);  % MODIFY

            %U1(:,i) = 0.0; % MODIFY
            
            %A1(i) = U1(1,i);
            %Q1(i) = U1(2,i);
            
            %========================
            % =============== FIXME
            % Lax-Wendroff with Strang splitting
            % i = 2, 3, ..., Nx-1
            
            % get 3-stencil data
            u_i = U(:,i);
            u_i_m1 = U(:,i-1); 
            u_i_p1 = U(:,i+1);
            
            % take two half steps
            
            u_mid_1 = (1/2)*(u_i + u_i_m1) + ...
                (dt/2)*( (flux(u_i_m1)-flux(u_i))/dx ...
                + ...
                (rhs(u_i)+rhs(u_i_m1))/2 ...
                );
            
            u_mid_2 = (1/2)*(u_i_p1 + u_i) + ...
                (dt/2)*( (flux(u_i)-flux(u_i_p1))/dx ...
                + ...
                (rhs(u_i_p1)+rhs(u_i))/2 ...
                );
            
            % advance one time step n -> n+1
            u1 = u_i - (dt/dx) * (flux(u_mid_2) - flux(u_mid_1)) + ...
                (dt/2)*(rhs(u_mid_2) + rhs(u_mid_1));
            % store solutions
            A1(i) = u1(1); 
            Q1(i) = u1(2);
            U1(:,i) = u1;
            
        end
    end
    A = A1; %update A, Q and U
    Q = Q1;
    U(1,:) = A;
    U(2,:) = Q;
 

  %plotting
   figure(1);
   width = 800; height = 600;
   set(figure(1),'Position',[15 15 width height]);
   if n==1 %For 1st layout
      subplot(2,2,1); plot(x,A);hold on;
      subplot(2,2,1); plot(x,ones(length(x))*R0*R0);hold off;
      axis([0 L  0.02 0.05]); title('Cross-sectional Area');
      xlabel('x (cm)'); ylabel('A (cm^2)');
      subplot(2,2,2); plot(x,Q); 
      axis([0 L -4 16]); title('Flow');
      xlabel('x (cm)'); ylabel('Q (cm^3/s)');
      subplot(2,2,3); plot(t,Inlet(t)/1333,'b');hold on; 
      subplot(2,2,3); plot(t,Outlet(t)/1333,'r');hold on;
      axis([0 T 50 200]); title('Proximal/Distal Pressure');
      ylabel('mmHg'); xlabel('t (sec)');
      subplot(2,2,4);plot(t(n),Q(1,N/2)/A(1,N/2)/3.14,'b-+');hold on;
      axis([0 T -20 160]); title('Velocity at mid-point');
      xlabel('t (sec)'); ylabel('V (cm/s)');
   end
   
   Nskip=200; %number of time steps to skip for plotting;
   if mod(n,Nskip)==1
     subplot(2,2,1); plot(x,A);hold on;
     subplot(2,2,1); plot(x,ones(length(x))*R0*R0);hold off;
     axis([0 L 0.02 0.05]); title('Cross-sectional Area');
      xlabel('x (cm)'); ylabel('A (cm^2)');
     subplot(2,2,2); plot(x,Q); 
     axis([0 L -4 16]); title('Flow');
     xlabel('x (cm)'); ylabel('Q (cm^3/s)');
     subplot(2,2,3);h2=plot([t(n) t(n)],[0, 250],'Color',[0 0 0]);
     subplot(2,2,4);plot(t(n),Q(1,N/2)/A(1,N/2)/3.14,'b-+');hold on;
     axis([0 T -20 160]); title('Velocity at mid-point');
      xlabel('t (sec)'); ylabel('V (cm/s)');
     drawnow;
     disp(t(n));
     
     %if n > (length(t)-1)/Ncycles
     figure(2);
     width = 800; height = 600;
     set(figure(2),'Position',[15 15 width height]);
     subplot(3,1,1);
     plot(t(n),Q(1,N/2)/A(1,N/2)/3.14,'*');hold on;
     axis([0.0 1.58 -20 160]); title('Velocity at mid-point');
     xlabel('t (sec)'); ylabel('V (cm/s)');
     subplot(3,1,2);
     plot(t(n),A(1,N/2),'*');hold on;
     axis([0.0 1.58 0.02 0.05]); title('Cross-sectional Area at mid-point');
     xlabel('t (sec)'); ylabel('A (cm^2)');
     
     subplot(3,1,3);
     plot(t(n),Q(1,N/2),'*');hold on;
     axis([0.0 2.5 2.5 15.0]); title('Flow Rate Q at mid-point');
     xlabel('t (sec)'); ylabel('Q (cm^3/s)');
     drawnow;
%     MovieFrame2(index2) = getframe(figure(2)) ;
%     index2=index2+1;
     %end
     
     
%Next 3 lines are generating movie frames     
     %MovieFrame1(index) = getframe(figure(1)) ;
     %set(h2,'Visible','off')
     %index=index+1;
%
   end
end

%Next 8 lines are saving the movie to an AVI file
%writerObj = VideoWriter('Case2.avi');
%writerObj.FrameRate = 15;
%open(writerObj);
%for i=1:length(MovieFrame1)
    %frame = MovieFrame1(i) ;    
    %writeVideo(writerObj, frame);
%end
%for i=1:length(MovieFrame2)
    %frame = MovieFrame2(i) ;    
    %writeVideo(writerObj, frame);
%end
%close(writerObj);
disp("maxflow reporting ...")
disp(['max flow= ', num2str(maxflow)])



%End of main procedure


%%%%%%%%%%%%%%%
%         Subroutines            %
%%%%%%%%%%%%%%%
% added more utilities ==========
% @author Hongli Zhao
function [l1,l2] = lambda_pair(u_star)
% a helper function that calculates the eigenvalues
% as directly given in the writeup
global alpha A0 G0 rho
%   input: 
%       u_star - 2x1 vector where we are linearizing
%   the system

%   output:
%       [l1,l2] - eigenvalues given directly with 
%   formulas in the writeup, for the system
aa = u_star(1);
qq = u_star(2);

l1 =  (alpha * qq)/aa + sqrt(alpha * (alpha-1) * (qq/aa)^2 ...
    + (G0/(2*rho))*(aa/A0)^(1/2));
l2 =  (alpha * qq)/aa - sqrt(alpha * (alpha-1) * (qq/aa)^2 ...
    + (G0/(2*rho))*(aa/A0)^(1/2));


end


function B = jacobi(u)
% a helper function that returns the jacobian matrix
% for the linearized system

global rho alpha

%   input: 
%       u - 2x1 vector, U_*, the solution around 
%   which we are linearizing the hyperbolic system
%
%
%   output:
%       B - 2x2 matrix representing the jacobian
%   of the flux

A_star = u(1); 
Q_star = u(2);
B = [0, 1;
    -alpha*(Q_star/A_star)^2 + (1/rho)*A_star*...
    p_prime(A_star), 2*alpha*Q_star/A_star
    ];
%global alpha G0 A0 rho
%Atru = u(1);
%Qtru = u(2);
%B = [0,1;-alpha*Qtru^2/Atru^2+sqrt(Atru/A0)*G0/(2*rho),...
%    2*alpha*Qtru/Atru];
end

function prime = p_prime(A)
% a helper function that calculates
% dp/dA
global A0 G0
prime = (G0/(2*sqrt(A0)))*(1/sqrt(A));
end
% ==============================



% calculating the flux term (in conservative form)
function f= flux(Utru) % F = [F1,F2]
global alpha rho R0 G0 A0 
Atru = Utru(1);
Qtru = Utru(2);
f = [Qtru; (alpha*Qtru^2/Atru)+(   G0/3/rho*((Atru/A0)^(3/2))*A0)];
end

%calculate the RHS
function s = rhs(Utru) 
global alpha nu 
Atru = Utru(1);
Qtru = Utru(2);
s = [0;(-2*nu*(alpha/(alpha-1))*(Qtru/Atru))];
end

% Return pressure boundary data at Inlet
function InP = Inlet(t)
global Xin T1
N = length(Xin);
if mod(N,2) == 0
    n = [0:(N/2 - 1) (-N/2):-1]';
else
    n = [0:floor(N/2) -floor(N/2):-1]';
end
InP = real(sum(Xin .* exp(2. * pi * 1i * n * mod(t,T1) / T1))) / N;
end

% Return pressure boundary data at Outlet
function OutP = Outlet(t)
global Xout T1
N = length(Xout);
if mod(N,2) == 0
    n = [0:(N/2 - 1) (-N/2):-1]';
else
    n = [0:floor(N/2) -floor(N/2):-1]';
end
OutP = real(sum(Xout .* exp(2. * pi * 1i * n * mod(t,T1) / T1))) / N;
end
