function [ci,ui,vi,pi,nui] = guess(n,dz,G,Mu,X0,Psi,cw,bb,m)
% this function assumes an initial guess for the concentration profile to
% allow enough concentration for the imposed sink
% n is the number of grid cells
% dz is the size of the grid cell
% Mu is the Munch number
% G is the ratio of hydrostatic pressure over osmotic pressure
% X0 is the ratio of xylem water potential scaling and osmotic pressure
% Psi is the nondimensional xylem water potential
% cw is the sucrose scaling in mol/m3
% bb is the coefficient that fits viscosity to concentration
% m is a case number, 0 == constant viscosity, 1 == variable viscosity

alpha = 6; % this factor is based on trial and error to allow enough concentration 
            % depending on the imposed sinks

% concentration
co         = alpha.*(G*dz.*( (1:n)' - 1/2 ) - X0.*Psi); % assumes concentration is balanced 
                                                        % by xylem water potential and hydrostatic pressure
switch m
    case 0
        nuo = ones(n,1);
    case 1
        nuo = exp( - bb(2)*cw.*(co - 1) ...
            - bb(3)*(cw^2).*(co.^2 - 1) ...
            - bb(4)*(cw^3).*(co.^3 - 1) ...
            - bb(5)*(cw^4).*(co.^4 - 1) );
end
[uo,vo,po] = VelocityS(co,nuo,Mu,Psi,X0,dz,n); % solves the pressure and velocity field
ci = co;
ui = uo;
vi = vo;
pi = po;
nui = nuo;



end