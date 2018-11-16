%% Defintion of model parameters

gam = 0.9;              % Number between 0 and 1
L = 10;                 % Length of a beam (in m)
m = 10;                 % Mass of a beam (in kg)
g = 9.81;               % Acceleration due to gravity (in kg/ms^2)
N = 10;                  % Number of times to average for each number of points

Points = 2:80;          % Range of numbers of points to try

Times = sparse(1,length(Points));

for i = 1:length(Points)
    disp('Number of points:')
    disp(Points(i))
    n = Points(i);

    x0 = ones(1,3*n);       % Starting vector  

    A_eq = sparse(5,3*n);  
    A_eq(1,1) = 1;          % Specifies linear constraint on x_0 
    A_eq(2,n+1) =  1;       % Specifies linear constraint on y_0 
    A_eq(3,2*n+1) = 1;      % Specifies linear constraint on z_0 
    A_eq(4,n) = 1;          % Specifies linear constraint on x_n 
    A_eq(5,2*n) = 1;        % Specifies linear constraint on y_n 

    b_eq = sparse(5,1);     % Specifies that all the variables above are zero (as constraints)
    b_eq(4) = n*gam*L;      % Specifies that x_n = gamma*L

    options = optimoptions('fmincon','Algorithm','interior-point');
    
%     options.MaxFunctionEvaluations = 1e6;
%     options.ConstraintTolerance = 1e-6;
%     options.StepTolerance = 1e-4;
%     options.MaxIterations = 1e6;
%     
    options.SpecifyObjectiveGradient = true;
    options.SpecifyConstraintGradient = true;

    
    TimeTrials = sparse(1,N);
    
    for j = 1:N
        tic;                    % Start clock to measure computational time
        [x, min, exitflag,output] = fmincon(@BeamGPE,x0,[],[],A_eq,b_eq,[],[],@BeamLength,options);
        min = m*g*min;          % Find overall minimum value (may not be needed)
        TimeTrials(j) = toc;    % Clock to end time, storing computation time in memory
    end
    Times(i) = mean(TimeTrials);
end

loglog(Points,Times,'r')

%% Function which returns the nonlinear constraints
% Returns c_eq, a vector containing the expressions in terms of lengths
% between successive points which 
function [c,c_eq,grad_c,grad_c_eq] = BeamLength(x)
L = 10;
gam = 0.9;
n = length(x)/3;

c_eq = sparse(1,n-1);
    for i = 1:n-1
        c_eq(i) = (x(i)-x(i+1))^2 + (x(i+n)-x(i+1+n))^2 + (x(i+2*n)-x(i+1+2*n))^2 - L^2;
    end
c = [];

grad_c = [];

grad_c_eq = sparse(n-1,3*n);
    for i = 1:n-1
        grad_c_eq(i,i) = 2*(x(i)-x(i+1));
        grad_c_eq(i,i+1) = 2*(x(i+1)-x(i));
        
        grad_c_eq(i,n+i) = 2*(x(n+i)-x(n+i+1));
        grad_c_eq(i,n+i+1) = 2*(x(n+i+1)-x(n+i));
        
        grad_c_eq(i,2*n+i) = 2*(x(2*n+i)-x(2*n+i+1));
        grad_c_eq(i,2*n+i+1) = 2*(x(2*n+i+1)-x(2*n+i));
    end
grad_c_eq = grad_c_eq';
end

%% Definition of function that returns the potential energy 
% Returns F, the gravitational potential energy (scalar)
% and gradF, the gradient of it (a vector)

function [F,gradF] = BeamGPE(x)

n = length(x)/3;
F = x;
F(1:n) = 0;
F(2*n+1:3*n) = 0;

F(n+1) = F(n+1)/2;
F(2*n) = F(2*n)/2;
F = sum(F);

gradF = zeros(1,length(x));
gradF(n+2:2*n-1) = 1;
gradF(n+1) = 0.5;
gradF(2*n) = 0.5;

end