function [Nu1,bb,rrs] = MLR_Viscosity(c,T,n)
%This function approximates the viscosity-sucrose relation using the multiple
%linear regression method

nu = Dynamic_Viscosity_T_Sucrose_CW(T,c);
Nu = log(nu);

%%% linear regression \nu = b2 + b1*c^2
% c2 = c.^2;
% b1 = ( mean(Nu)*mean(c2.*c2) - mean(c2.*Nu')*mean(c2) )/(mean(c2.*c2) - mean(c2)*mean(c2));
% b2 = (mean(Nu) - b1)/mean(c2);
% bb = [b2 0 b1 0];
% Nu1 = c2.*b2 + b1;
% rrs = sqrt(sum((Nu - Nu1').^2)/n);

%%% Multiple linear regression \nu = b1 + b2*c + b3*c^2 + b4*c^3
c2 = c.^2;
c3 = c.^3;
c4 = c.^4;
xx  = [ones(n,1),c,c2,c3,c4];
mm  = inv(mtimes(xx',xx));
mm1 = mtimes(mm,xx');
bb  = mtimes(mm1,Nu');
Nu1 = bb(1) + bb(2).*c + bb(3).*c2 + bb(4).*c3 + bb(5).*c4;
rrs = sqrt(sum((Nu - Nu1').^2)/n);

end