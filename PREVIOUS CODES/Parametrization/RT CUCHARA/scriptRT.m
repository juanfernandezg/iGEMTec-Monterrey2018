%% Parametros inicial
r0 = rand(1);
C0 = rand(1)*.9;
d0 = r0;
b0 = rand(1)*.7;
x0 = [r0,C0,d0,b0];
%x0=x;
%x0 = [1E-15,1E-15,1E-15,1E-15];
format long
prec = 1e-45;
options = optimoptions('lsqnonlin','OptimalityTolerance',prec,'FunctionTolerance',prec,'StepTolerance',prec);
load('data2.mat','data2')
data=data2;
options.MaxFunctionEvaluations = 5000;
options.MaxIterations = 5000;

fun = @(s)errores(s)-data;

x = lsqnonlin(fun,x0,[0,0,0,0],[],options)
res=errores(x);
es=sum(sum((data-res).^2))
