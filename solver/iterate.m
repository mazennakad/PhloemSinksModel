function [ck,uk,vk,pk,nuk,iteration] = iterate(n,dz,Mu,X0,Pe,Sl,Ss,Sr,es,...
        Psi,co,uo,vo,po,nuo,cw,bb,k,k1,s,dt,cp,Ss1,Ss2)
% This function performs the newtons iteration with subscript o denoting
% previous and subscript k denoting current
check_C = true; % for stopping the iteration
iteration = 1;  % for the number of iterations
resk = 1;

while check_C || resk > 1e-7
    % Set the new variables
    if iteration ~= 1
        co = ck;
        uo = uk;
        vo  = vk;
        po  = pk;
        nuo = nuk;
    end

    [vector_S,Sk,resk] = newton(n,dz,Mu,X0,Pe,Sl,Ss,Sr,es,...
        Psi,co,uo,vo,po,nuo,cw,bb,k,k1,s,dt,cp,Ss1,Ss2);
    check_C = check_convergence(vector_S,Sk,length(Sk));
    pk  = Sk(1:n);
    uk  = Sk(n+1:2*n-1);
    vk  = Sk(2*n:3*n-1);
    ck  = Sk(3*n:4*n-1);
    nuk = Sk(4*n:end);

    if iteration > 100000
        break;
    end
    iteration = iteration + 1; % to set a maximum number of iterations
    
end

% fprintf('iteration done: %d \n',iteration-1)


end
