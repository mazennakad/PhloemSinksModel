function [VP,VU,VV,VC,VN,F3] = BuildV(n,dz,v,p,nu)
% This function builds the Jacobian block matrix for the V-equation
% n is the number of grid nodes
% dz is the grid size
% v is the nondimensional radial velocity, nx1
% p is the nondimensional dynamic pressure, nx1
% nu is the nondimensional inverse of viscosity, nx1

% VP is the Jacobian block matrix of size nxn for the P-variable
% VU is the Jacobian block matrix of size nx(n-1) for the U-variable
% VV is the Jacobian block matrix of size nxn for the V-variable
% VC is the Jacobian block matrix of size nxn for the C-variable
% VN is the Jacobian block matrix of size nxn for the Nu-variable
% F3 is the equality (residual) for the V-equation

% P variable
diagP = - nu(2:n-1)./(12*(dz^2));
uppP  = (nu(3:n) - nu(1:n-2))./(96*(dz^2)) ...
      + nu(2:n-1)./(24*(dz^2));
lowP  = - (nu(3:n) - nu(1:n-2))./(96*(dz^2)) ...
      + nu(2:n-1)./(24*(dz^2));
diagP = [0; diagP; 0];
uppP  = [0; uppP];
lowP  = [lowP; 0];
VP    = diag(diagP) + diag(lowP,-1) + diag(uppP,1);

% U variable
VU    = zeros(n,n-1); 
% V variable 
VV    = diag(ones(1,n)); 
% C variable
VC    = zeros(n);

% Nu variable
diagN = (p(3:n) - 2.*p(2:n-1) + p(1:n-2))./(24*(dz^2));
uppN  = (p(3:n) - p(1:n-2))./(96*(dz^2));
lowN  = - uppN;
diagN = [0; diagN; 0];
uppN  = [0; uppN];
lowN  = [lowN; 0];
VN    = diag(diagN) + diag(lowN,-1) + diag(uppN,1);

% Equality
F3 = v(2:n-1) ...
    + (nu(3:n) - nu(1:n-2)).*(p(3:n) - p(1:n-2))./(96*(dz^2))...
    + nu(2:n-1).*(p(3:n) - 2.*p(2:n-1) + p(1:n-2))./(24*(dz^2));
F3 = [0; F3; 0];

end