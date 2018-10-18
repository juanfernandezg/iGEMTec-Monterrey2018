function [difs] = errores(x)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


%Poblaci�n inicial
p1 = 2576249518; % 0mM
p2 = 2696489168; % 0.1mM
p3 = 2275650393; % 0.5mM
p4 = 2290162075; % 1.0mM

prot0=[0.690397350993;
0.706953642384;
0.392384105960;
0.442052980132];
prot1=prot0(1);
prot2=prot0(2);
prot3=prot0(3);
prot4=prot0(4);

r=x(1);
C=x(2);
d=x(3);
b=x(4);

%Condiciones iniciales
I1 = 0;
I2 = 0.1;
I3 = 0.5;
I4 = 1;

format long
h = .05;       %step
tlim = 6;  % horas

dydt = @(t,y) [8.63383E+08; %Crecimiento poblacional [mL^-1]
               ((r.*I1)./(C + I1) + b).*1 - d.*y(2); %Crecimiento prot�ico
               1.29998E+09; 
               ((r.*I2)./(C + I2) + b).*1 - d.*y(4); 
               4.87913E+07.*t*2 + 7.83194E+08;  
               ((r.*I3)./(C + I3) + b).*1 - d.*y(6); 
               1.62679E+08.*t*2 + 4.24091E+08; 
               ((r.*I4)./(C + I4) + b).*1 - d.*y(8)
               ];

[t, y] = ode45(dydt, 1:h:tlim, [p1;prot1;p2;prot2;p3;prot3;p4;prot4]);
difs = [y(find(t==1),[2,4,6,8]);y(find(t==2),[2,4,6,8]);y(find(t==4),[2,4,6,8]);y(find(t==6),[2,4,6,8])]';

end

