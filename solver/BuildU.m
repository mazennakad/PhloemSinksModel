function [UP,UU,UV,UC,UN,F2] = BuildU(n,dz,u,p,nu)
% This function builds the Jacobian block matrix for the U-equation
% n is the number of grid nodes
% dz is the grid size
% u is the nondimensional axial velocity, (n-1)x1
% p is the nondimensional dynamic pressure, nx1
% nu is the nondimensional inverse of viscosity, nx1

% UP is the Jacobian block matrix of size (n-1)xn for the P-variable
% UU is the Jacobian block matrix of size (n-1)x(n-1) for the U-variable
% UV is the Jacobian block matrix of size (n-1)xn for the V-variable
% UC is the Jacobian block matrix of size (n-1)xn for the C-variable
% UN is the Jacobian block matrix of size (n-1)xn for the Nu-variable
% F2 is the equality (residual) for the U-equation

% P variable
UP = diag( - (nu(2:n) + nu(1:n-1))./(24*dz) ) ... 
    + diag( (nu(2:n-1) + nu(1:n-2))./(24*dz),1 );
UP = [UP zeros(n-1,1)];
UP(n-1,n) = (nu(n) + nu(n-1))/(24*dz);

% U variable
UU = diag(ones(1,n-1)); 
% V variable 
UV = zeros(n-1,n);
% C variable
UC = zeros(n-1,n);
% Nu variable
UN = diag( (p(2:n) - p(1:n-1))./(24*dz) ) ...
    + diag( (p(2:n-1) - p(1:n-2))./(24*dz),1 );
UN = [UN zeros(n-1,1)];
UN(n-1,n) = (p(n) - p(n-1))/(24*dz);

% Equality
F2 = u + (nu(2:n) + nu(1:n-1)).*(p(2:n) - p(1:n-1))./(24*dz);
end