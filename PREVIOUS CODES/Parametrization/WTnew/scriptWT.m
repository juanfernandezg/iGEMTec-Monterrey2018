%% Parametros inicial
r0 = rand(1);
C0 = rand(1);
d0 = rand(1);
b0 = rand(1);
%x0 = [r0,C0,d0,b0];
x0 = x;
%x0 = [1E-15,1E-15,1E-15,1E-15];
format long
prec = 1e-20;
options = optimoptions('lsqnonlin','OptimalityTolerance',prec,'FunctionTolerance',prec,'StepTolerance',prec);
load('data2.mat','data')
options.MaxFunctionEvaluations = 400;
options.MaxIterations = 400;
fun = @(s)errores(s)-data;

x = lsqnonlin(fun,x0,[],[],options)

