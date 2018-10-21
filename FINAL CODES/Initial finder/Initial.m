%% Initial parameters
x0 = [0,0,.5,0,0,0,0,0,0,0,0,0];
%x0=x;
%x0 = [1E-15,1E-15,1E-15,1E-15];

%% ODE system sum of square errors (sse) minimization
% Minimization configuration
format long
prec = 1e-30;
options = optimoptions('lsqnonlin','OptimalityTolerance',prec,'FunctionTolerance',prec,'StepTolerance',prec);
options.MaxFunctionEvaluations = 4000;
options.MaxIterations = 4000;

% Experimental data
load('final.mat','final')

% Error function
error = @(s)aprox(s)-final;

% Errors minimization
x = lsqnonlin(error,x0,[],[],options)
% Differential Equations Systems Solution
res=aprox(x);
%Square errors
se=(final-res).^2;
%Sum of square errors
sse=sum(sum(se))


% Results
a0 = x(1);
a1 = x(2);
a2 = x(3);
a3 = x(4);
a4 = x(5);
a5 = x(6);
a6 = x(7);
a7 = x(8);
a8 = x(9);
a9 = x(10);
a10 = x(11);
a11 = x(12);

yapprox = a0 + a1.*t + a2.*t.^2 + a3.*t.^3 + a4.*t.^4 + a5.*t.^5;% + 6*a6.*t.^5 + 7*a7.*t.^6 + 8*a8.*t.^7 + 9*a9.*t.^8 + 10*a10.*t.^9 + 11*a11.*t.^10;%Taylor
figure(23)
plot(t,yapprox)
