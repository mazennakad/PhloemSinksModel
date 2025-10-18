function [c,u,nu,v,p] = initial(n,dz,Mu,G,X0,Pe,Sl,Ss,Sr,es,Psi,...
            cw,bb,m,k1,Ss1,Ss2,dt)
% This function assumes an initial profile for the variables
% guess function assumes an initial profile for concentration and solves
% for other variables
% iterate function performs newton's method for the steady state
% conditions, i.e. no time derivative

% [c,u,v,p,nu] = guess(n,dz,G,Mu,X0,Psi,cw,bb,m);
[co,uo,vo,po,nuo] = guess(n,dz,G,Mu,X0,Psi,cw,bb,m);
[c,u,v,p,nu,~]    = iterate(n,dz,Mu,X0,Pe,Sl,Ss,Sr,es,...
                            Psi,co,uo,vo,po,nuo,cw,bb,m,k1,0,dt,zeros(n,1),Ss1,Ss2);

end