function [check_S] = check_convergence(vector_S,Sk,N)
% This function checks if the iteration has converged based on root mean
% square

rms = sqrt(sum((vector_S - Sk).^2)/N);
if rms <= 1e-3
    check_S = false;
else
    check_S = true;
end
end

