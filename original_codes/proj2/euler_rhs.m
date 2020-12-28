function [fr,fru,frv,frE] = euler_rhs(r, ru, rv, rE, h)
% computes the right-hand side of the discretized divergence
% there is a negative sign

% inputs:
%   - r, ru, rv, rE: input (N+1)x(N+1) matrices from which we 
% would like to calculate the divergence for
%   - h: step size of your discretization, to be fed in when 
% approximating divergence

% outputs: 
%   - fr, fru, fry, frE: the four divergence components

    % Get all 8 fluxes:
    [Frx,Fry,Frux,Fruy,Frvx,Frvy,FrEx,FrEy] = euler_fluxes(r, ...
    ru, rv, rE);

    fr = -compact_div(Frx,Fry,h);
    fru = -compact_div(Frux,Fruy,h);
    frv = -compact_div(Frvx,Frvy,h);
    frE = -compact_div(FrEx,FrEy,h);

end