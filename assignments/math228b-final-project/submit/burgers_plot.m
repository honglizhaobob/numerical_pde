% Generates Plots for various cases in problem 1
h = 0.01;
x = -1.5:h:1.5;
t = 0:h:0.5;
Xi = -1:0.04:1;

mode = 6; % 1, 2, 3, 4, 5, 6
if mode == 1
   % t=0, plot u
   figure(1)
   plot(x, sin(x*pi), 'LineWidth', 1.2, 'Color', 'r');
   title('plot of solution for t=0')
   xlabel('x'); ylabel('u(x,0)');
   grid on
end

if mode == 2
    % characteristic x-t curve, envelope
   figure(2)
   grid on
   %xt = sin(acos(-1./(pi*t))) + (1/pi)*(acos(-1./(pi*t)));
   %plot(t,xt)
   % plot a sample of curves by varying xi
   for i = 1:length(Xi)
       hold on
       xi = Xi(i);
      %xt = sin(pi*xi).*t + xi;
       tx = (x - xi)/sin(pi*xi);
       plot(x,tx,'LineWidth',1.2)
   end
   title('char lines x-t');
   xlim([-1.5 1.5]); ylim([0 0.5]);
   xlabel('x'); ylabel('t(x)');
end

if mode == 3
    % implicit solution t=0.25
    figure(1)
    grid on
    f = @(x,u) u - sin(pi*(x-0.25*u));
    fimplicit(f,[-1,1],'LineWidth',1.2,'Color','r')
    title('plot of solution for t=0.25')
    xlabel('x'); ylabel('u(x,0.25)');
   grid on
    
end

if mode == 4
    % implicit solution t = 0.5
    figure(1)
    grid on
    f = @(x,u) u - sin(pi*(x-0.5*u));
    fimplicit(f,[-1.5 1.5],'LineWidth',1.2,'Color','r')
    title('plot of solution for t=0.5')
    xlabel('x'); ylabel('u(x,0.5)');
    grid on
end

if mode == 5
    % plotting basis functions
    figure(1)

    phi_10 = -2*x-1;
    phi_11 = 2*x+2;
    
    phi_20 = -2*x;
    phi_21 = 2*x+1;
    
    phi_30 = -2*x+1;
    phi_31 = 2*x;
    
    phi_40 = -2*x+2;
    phi_41 = 2*x-1;
    plot(x,phi_10,x,phi_11,x,phi_20,x,phi_21, ...
        x,phi_30,x,phi_31,x,phi_40,x,phi_41,'LineWidth',1.2)
    ylim([0 1])
    grid on
end

if mode == 6
   % symbolic computations for matrices
    syms x
    phi_10 = -2*x-1;
    phi_11 = 2*x+2;

    phi_20 = -2*x;
    phi_21 = 2*x+1;

    phi_30 = -2*x+1;
    phi_31 = 2*x;

    phi_40 = -2*x+2;
    phi_41 = 2*x-1;
    
    k1 = [-1,-1/2]; k2 = [-1/2,0]; k3 = [0,1/2]; k4 = [1/2,1];
    
    M1 = [int(phi_10^2, k1), int(phi_10*phi_11, k1);
          int(phi_11*phi_10, k1), int(phi_11^2, k1)];
    M2 = [int(phi_20^2, k2), int(phi_20*phi_21, k2);
          int(phi_21*phi_20, k2), int(phi_21^2, k2)];
    M3 = [int(phi_30^2, k3), int(phi_30*phi_31, k3);
          int(phi_31*phi_30, k3), int(phi_31^2, k3)];
    M4 = [int(phi_40^2, k4), int(phi_40*phi_41, k4);
          int(phi_41*phi_40, k4), int(phi_41^2, k4)];
    C1 = [int(phi_10*diff(phi_10),k1), int(phi_10*diff(phi_11),k1);
        int(phi_11*diff(phi_10),k1), int(phi_11*diff(phi_11),k1)];
    C2 = [int(phi_20*diff(phi_20),k2), int(phi_20*diff(phi_21),k2);
        int(phi_21*diff(phi_20),k2), int(phi_21*diff(phi_21),k2)];
    C3 = [int(phi_30*diff(phi_30),k3), int(phi_30*diff(phi_31),k3);
        int(phi_31*diff(phi_30),k3), int(phi_31*diff(phi_31),k3)];
    C4 = [int(phi_40*diff(phi_40),k4), int(phi_40*diff(phi_41),k4);
        int(phi_41*diff(phi_40),k4), int(phi_41*diff(phi_41),k4)];
     
end