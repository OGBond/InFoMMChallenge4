%% Function which returns the nonlinear constraints
% Returns c_eq, a vector containing nonlinear expressions which are to 
% be equal to zero (i.e. length between consecutive endpoints of beams
% minus length squared.

function [c,c_eq,grad_c,grad_c_eq] = BeamLength(x)

L = 1;
n = length(x)/2;

c_eq = zeros(1,n-1);
    for i = 1:n-1
        c_eq(i) = (x(i)-x(i+1))^2 + (x(i+n)-x(i+1+n))^2 - L^2;
    end
c = [];

grad_c = [];

grad_c_eq = zeros(n-1,2*n);
    for i = 1:n-1
        grad_c_eq(i,i) = 2*(x(i)-x(i+1));
        grad_c_eq(i,i+1) = 2*(x(i+1)-x(i));
        
        grad_c_eq(i,n+i) = 2*(x(n+i)-x(n+i+1));
        grad_c_eq(i,n+i+1) = 2*(x(n+i+1)-x(n+i));
        
    end
grad_c_eq = grad_c_eq';
end
