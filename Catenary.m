%% Defintion of model parameters

L = 1;
m = 10;
g = 9.81;
gam = 0.5;              % Number between 0 and 1;
N = 3;                  % Number of times to average for each number of points

%% Declaration of plot features

Methods = ["active-set";"interior-point";"sqp";"sqp-legacy"];
Colors = ["r";"g";"b";"m"];
Colors = Colors + "x";

%% Storage space for key information

Points = 150:5:200;          % Range of numbers of points to try
Times = zeros(length(Methods),length(Points));

%% Solves the optimisation problem using each method [

for k = 1:size(Methods)
    
    disp(Methods(k))

    for i = 1:length(Points)
        disp('Number of points:')
        disp(Points(i))
        n = Points(i);

        x0 = rand(1,2*n);       % Starting vector, with randomly generated
                                % between 0 and 1
        A_eq = zeros(4,2*n);  
        A_eq(1,1) = 1;          % Specifies linear constraint on x_0 
        A_eq(2,n+1) =  1;       % Specifies linear constraint on y_0 
        A_eq(3,n) = 1;          % Specifies linear constraint on x_n 
        A_eq(4,2*n) = 1;        % Specifies linear constraint on y_n 

        b_eq = zeros(4,1);     % Specifies that all the variables above are zero (as constraints)
        b_eq(3) = n*gam*L;      % Specifies that x_n = gamma*L

        options = optimoptions('fmincon');
        
        options.Algorithm = Methods(k);
        options.MaxFunctionEvaluations = 1e6;
        options.ConstraintTolerance = 1e-6;
        options.StepTolerance = 1e-6;
        options.MaxIterations = 1e9;
        options.FunctionTolerance = 1e-9;

        options.SpecifyObjectiveGradient = true;
        options.SpecifyConstraintGradient = true;

        TimeTrials = zeros(1,N);

        for j = 1:N
            tic;                    % Start clock to measure computational time
            [x, minimum, exitflag,output] = fmincon(@BeamGPE,x0,[],[],A_eq,b_eq,[],[],@BeamLength,options);
            minimum = m*g*minimum;          % Find overall minimum value (may not be needed)
            TimeTrials(j) = toc;    % Clock to end time, storing computation time in memory
        end
        
        if exitflag > 0
            Times(k,i) = mean(TimeTrials);
        else
            Times(k,i) = NaN;
        end
    end
end

%% Plot log-log graphs of computation times against points for each method

for k = 1:length(Methods)
    plot(Points,Times(k,:),Colors(k))
    hold on
end
set(gca, 'XScale', 'log', 'YScale', 'log');
legend(Methods)
hold off

%% Calculate approximate polynomial complexities from these points
ActiveSetP = polyfit(log(Points),log(Times(1,:)),1);
ActiveSetComplexity = ActiveSetP(1)

InteriorPointP = polyfit(log(Points),log(Times(2,:)),1);
InteriorPointComplexity = InteriorPointP(1)

SQP_P = polyfit(log(Points),log(Times(3,:)),1);
SQP_Complexity = SQP_P(1)

SQPLegacyP = polyfit(log(Points),LogTimes(4,:),1);
SQPLegacyComplexity = SQPLegacyP(1)

%% Plot solution

% plot(x(1:n),x(n+1:2*n))
% % xlim([min(x(1:n)),max(x(1:n))])
% xlabel('$x$','Interpreter','LaTeX','FontSize',15)
% ylabel('$y$','Interpreter','LaTeX','FontSize',15)
