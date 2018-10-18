function dydt = DESystemWT2(t,y)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%Condiciones iniciales
I1 = 0.1;
I2 = 0.1;
I3 = 0.1;
I4 = 0.1;

%Parametros a encontrar
r = .00001;
C = 0.0000001;
d = 0.0000001;

dydt = [8.89065E+08; %Crecimiento poblacional [mL^-1]
        (r.*I1./y(1))./(C + I./y(1)) - d*y(2); %Crecimiento protéico
        1.299984*10^9; %Crecimiento poblacional [mL^-1]
        (r.*I2./y(3))./(C + I./y(3)) - d*y(4); %Crecimiento protéico
        (r.*I3./y(5))./(C + I./y(5)) - d*y(6); %Crecimiento protéico
        (r.*I4./y(7))./(C + I./y(7)) - d*y(8); %Crecimiento protéico
        ];
end

