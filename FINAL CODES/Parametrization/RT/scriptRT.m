%% Initial parameters
r0 = rand(1)*.10;
C0 = rand(1)*.09;
d0 = r0;
b0 = rand(1)*.07;
n0=1;
x0 = [r0,C0,d0,b0];
%x0=x;
%x0 = [1E-15,1E-15,1E-15,1E-15];

%% ODE system sum of square errors (sse) minimization
% Minimization configuration
format long
prec = 1e-20;
options = optimoptions('lsqnonlin','OptimalityTolerance',prec,'FunctionTolerance',prec,'StepTolerance',prec);
options.MaxFunctionEvaluations = 400;
options.MaxIterations = 400;

% Experimental data
load('data.mat','data')

% Error function
error = @(s)aprox(s)-data;

% Errors minimization
x = lsqnonlin(error,x0,[0,0,0,0],[1000,1,1000,1000],options)
% Differential Equations Systems Solution
res=aprox(x);
%Square errors
se=(data-res).^2
%Sum of square errors
sse=sum(sum(se))
