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



