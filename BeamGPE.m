%% Definition of function that returns the potential energy 
% Returns F, the gravitational potential energy (scalar)
% and gradF, the gradient of it (a vector)

function [F,gradF] = BeamGPE(x)

m = 10;
g = 9.81;
n = length(x)/2;
F = x;
F(1:n) = 0;

F(n+1) = F(n+1)/2;
F(2*n) = F(2*n)/2;
F = sum(F);

gradF = zeros(1,length(x));
gradF(n+2:2*n-1) = 1;
gradF(n+1) = 0.5;
gradF(2*n) = 0.5;

end