function [ui,vi,pi] = VelocityS(ci,nui,Mu,Psi,X0,dz,n)
% This function solves the initial profile for veloctiy and pressure
% ui is the axial velcoity
% vi is the radial velocity
% pi is the dynamic pressure
% ci is the concentration
% nui is the inverse of viscosity
% Mu is the Munch number
% Psi is the xylem water potential
% X0 is the scaling of xylem water potential over osmosis
% dz is the grid size
% n is the number of grid cells


% pressure
diagP      = - nui(2:n-1).*2/(12*(dz^2)) - Mu;
diagP1     = - Mu; 
diagP2     = - Mu;
diagP      = [diagP1; diagP; diagP2];
uppP       = (nui(3:n) - nui(1:n-2))./(48*(dz^2)) + nui(2:n-1)./(12*(dz^2));
uppP       = [0; uppP;0];
lowP       = - (nui(3:n) - nui(1:n-2))./(48*(dz^2)) + nui(2:n-1)./(12*(dz^2));
lowP       = [0;lowP; 0];

eqP        = - ci(2:n-1) - X0.*Psi(2:n-1);
eqP1       = - ci(1)  -  X0*Psi(1);
eqP2       = - ci(n) - X0*Psi(n);
eqP        = [eqP1;eqP;eqP2];
pi         = Thomas(lowP,diagP,uppP,eqP);
pi         = pi';

% axial velocity
ui         = (nui(2:n) + nui(1:n-1)).*( pi(1:n-1) - pi(2:n) )./(24*dz);

% radial velocity
vi         = zeros(n,1);
vi(2:n-1)  = (1/2).*(ui(2:n-1) - ui(1:n-2))./dz;

end