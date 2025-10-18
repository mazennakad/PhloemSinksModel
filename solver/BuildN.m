function [NP,NU,NV,NC,NN,F5] = BuildN(n,cw,bb,c,nu)
% This function builds the Jacobian block matrix for the Nu-equation
% n is the number of grid nodes
% cw is the sucrose scaling in mol/m3
% bb is 5x1 vector that has the coefficient of multiple linear regression
% for viscosity as a function of concentration
% c is the nondimensional concentration, nx1
% nu is the nondimensional inverse of viscosity, nx1

% NP is the Jacobian block matrix of size nxn for the P-variable
% NU is the Jacobian block matrix of size nx(n-1) for the U-variable
% NV is the Jacobian block matrix of size nxn for the V-variable
% NC is the Jacobian block matrix of size nxn for the C-variable
% NN is the Jacobian block matrix of size nxn for the Nu-variable
% F5 is the equality (residual) for the Nu-equation

% P variable
NP = zeros(n);
% U variable
NU = zeros(n,n-1);
% V variable
NV = zeros(n);
% C variable 
n1 =  - ( - bb(2)*cw - 2*bb(3)*(cw^2).*c ...
    - 3*bb(4)*(cw^3).*c.^2 - 4*bb(5)*(cw^4).*c.^3 ) ...
   .*exp( - bb(2)*cw.*(c - 1) - bb(3)*(cw^2).*(c.^2 - 1) ...
    - bb(4)*(cw^3).*(c.^3 - 1) - bb(5)*(cw^4).*(c.^4 - 1) ) ;
NC = diag(n1);
% Nu variable
NN = diag(ones(1,n));

% Equality
F5 = nu - exp( - bb(2)*cw.*(c - 1) - bb(3)*(cw^2).*(c.^2 - 1) ...
    - bb(4)*(cw^3).*(c.^3 - 1) - bb(5)*(cw^4).*(c.^4 - 1) );
end