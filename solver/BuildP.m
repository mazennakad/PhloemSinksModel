function [PP,PU,PV,PC,PN,F1] = BuildP(n,dz,Mu,X0,Psi,c,p,nu)
% This function builds the Jacobian block matrix for the P-equation
% n is the number of grid nodes
% dz is the grid size
% Mu is the Munch number
% X0 is the ratio of xylem water potential scale (leaf level) over the
% osmotic potential scale
% Psi is the nondimensional xylem water potential, nx1
% c is the nondimensional concentration, nx1
% p is the nondimensional dynamic pressure, nx1
% nu is the nondimensional inverse of viscosity, nx1

% PP is the Jacobian block matrix of size nxn for the P-variable
% PU is the Jacobian block matrix of size nx(n-1) for the U-variable
% PV is the Jacobian block matrix of size nxn for the V-variable
% PC is the Jacobian block matrix of size nxn for the C-variable
% PN is the Jacobian block matrix of size nxn for the Nu-variable
% F1 is the equality (residual) for the P-equation

% P variable
diagP  = - nu(2:n-1).*2/(12*(dz^2)) - Mu;
diagP1 = - Mu; 
diagP2 = - Mu;
diagP  = [diagP1; diagP; diagP2];
uppP   = (nu(3:n) - nu(1:n-2))./(48*(dz^2)) + nu(2:n-1)./(12*(dz^2));
uppP   = [0; uppP];
lowP   = - (nu(3:n) - nu(1:n-2))./(48*(dz^2)) + nu(2:n-1)./(12*(dz^2));
lowP   = [lowP; 0];
PP     = diag(diagP) + diag(uppP,1) + diag(lowP,-1);

% U variable
PU     = zeros(n,n-1);
% V variable
PV     = zeros(n);
% C variable
PC     = diag(ones(1,n));

% Nu variable
diagN  = ( p(3:n) - 2.*p(2:n-1) + p(1:n-2) )./(12*(dz^2));
diagN  = [0; diagN; 0];
uppN   = (p(3:n) - p(1:n-2))./(48*(dz^2));
lowN   = - uppN;
uppN   = [0; uppN];
lowN   = [lowN; 0];
PN     = diag(diagN) + diag(uppN,1) + diag(lowN,-1);

% Equality
F1     = (nu(3:n) - nu(1:n-2)).*(p(3:n) - p(1:n-2))./(48*(dz^2)) ...
       + nu(2:n-1).*( p(3:n) - 2.*p(2:n-1) + p(1:n-2) )./(12*(dz^2)) ...
       - Mu.*p(2:n-1) + c(2:n-1) + X0.*Psi(2:n-1);
F11    = - Mu*p(1) + c(1) + X0*Psi(1);
F12    = - Mu*p(n) + c(n) + X0*Psi(n);
F1     = [F11; F1; F12];


end