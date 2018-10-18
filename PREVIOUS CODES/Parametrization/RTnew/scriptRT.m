%% Parametros inicial
r0 = rand(1)*2e-5;
C0 = rand(1)*1e-6;
d0 = rand(1);
b0 = rand(1)*1e-10;
%x0 = [r0,C0,d0,b0];
x0=x;
%x0 = [1E-15,1E-15,1E-15,1E-15];
format long
prec = 1e-45;
options = optimoptions('lsqnonlin','OptimalityTolerance',prec,'FunctionTolerance',prec,'StepTolerance',prec);
load('data.mat','data')
options.MaxFunctionEvaluations = 5000;
options.MaxIterations = 5000;

fun = @(s)errores(s)-data;

x = lsqnonlin(fun,x0,[0,0,0,0],[],options)

