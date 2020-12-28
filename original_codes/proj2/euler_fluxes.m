function [Frx,Fry,Frux,Fruy,Frvx,Frvy,FrEx,FrEy] = euler_fluxes(r, ...
    ru, rv, rE)
% inputs: r, ru, rv, rE (vector components of u), each is a scalar
% (written in general case so easily switched if needed)
% outputs: flux functions of F, scalar, at (x,y)
    % on phyical domain
    u = ru./r;
    v = rv./r;
    p = ((7/5)-1)*(rE - (ru.*u + rv.*v)/2);   

    Frx = ru; 
    Fry = rv;
    Frux = ru.*u + p;
    Fruy = ru.*v;
    Frvx = Fruy;
    Frvy = rv.*v + p;
    FrEx = u.*(rE + p);
    FrEy = v.*(rE + p);
end