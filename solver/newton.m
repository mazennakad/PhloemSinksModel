function [vector_S,Sk,resk] = newton(n,dz,Mu,X0,Pe,Sl,Ss,Sr,es,Psi,...
    c,u,v,p,nu,cw,bb,k,k1,s,dt,cp,Ss1,Ss2)
% This function solves Newton's method
% m is the Jacobian matrix
% scaling is the regularization parameter
% sol is the update vector

[PP,PU,PV,PC,PN,F1] = BuildP(n,dz,Mu,X0,Psi,c,p,nu);
[UP,UU,UV,UC,UN,F2] = BuildU(n,dz,u,p,nu);
[VP,VU,VV,VC,VN,F3] = BuildV(n,dz,v,p,nu);
[CP,CU,CV,CC,CN,F4] = BuildC(n,dz,Pe,Sl,Ss,Sr,es,c,u,v,s,dt,cp,k1,Ss1,Ss2);

switch k

    case 0
        NP = zeros(n);
        NU = zeros(n,n-1);
        NV = zeros(n);
        NC = zeros(n);
        NN = diag(ones(1,n));
        F5 = zeros(n,1);
    case 1
        [NP,NU,NV,NC,NN,F5] = BuildN(n,cw,bb,c,nu);    
end

m        = sparse([sparse(PP) sparse(PU) sparse(PV) sparse(PC) sparse(PN);...
            sparse(UP) sparse(UU) sparse(UV) sparse(UC) sparse(UN);...
            sparse(VP) sparse(VU) sparse(VV) sparse(VC) sparse(VN);...
            sparse(CP) sparse(CU) sparse(CV) sparse(CC) sparse(CN);...
            sparse(NP) sparse(NU) sparse(NV) sparse(NC) sparse(NN)]);
F        = - [F1;F2;F3;F4;F5];

scaling = 1e-6;
sol = (m+scaling.*(speye(length(F))))\F;

vector_S = [p; u; v; c; nu];
Sk       = vector_S + sol;

% Check the residual
pk  = Sk(1:n);
uk  = Sk(n+1:2*n-1);
vk  = Sk(2*n:3*n-1);
ck  = Sk(3*n:4*n-1);
nuk = Sk(4*n:end);

[~,~,~,~,~,F1] = BuildP(n,dz,Mu,X0,Psi,ck,pk,nuk);
[~,~,~,~,~,F2] = BuildU(n,dz,uk,pk,nuk);
[~,~,~,~,~,F3] = BuildV(n,dz,vk,pk,nuk);
[~,~,~,~,~,F4] = BuildC(n,dz,Pe,Sl,Ss,Sr,es,ck,uk,vk,s,dt,cp,k1,Ss1,Ss2);
switch k

    case 0
        F5 = zeros(n,1);
    case 1
        [~,~,~,~,~,F5] = BuildN(n,cw,bb,c,nu);    
end

Fk        = -[F1;F2;F3;F4;F5];
resk = norm(Fk,inf);
end