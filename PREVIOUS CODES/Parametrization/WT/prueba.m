
%Población inicial
p1 = 1.72437E+09; % 0mM
p2 = 1.41015E+09; % 0.1mM
p3 = 1.41105E+09; % 0.5mM
p4 = 1.57027E+09; % 1.0mM
prot0 = 4.062415619E-10; 
r=x(1);
C=x(2);
d=x(3);
b=x(4);

%Condiciones iniciales
I1 = 0;
I2 = 0.1;
I3 = 0.5;
I4 = 1;

% data = [4.06000000000000e-10,2.67985000000000e-10,1.76659000000000e-10,1.21718000000000e-10,1.68388000000000e-10;4.06000000000000e-10,2.62176000000000e-10,1.33167000000000e-10,1.16143000000000e-10,1.59556000000000e-10;4.06000000000000e-10,1.72427000000000e-10,2.13603000000000e-10,1.20391000000000e-10,1.29818000000000e-10;4.06000000000000e-10,1.93023000000000e-10,2.22328000000000e-10,1.71893000000000e-10,1.01983000000000e-10];

format long
h = .001;       %step
tlim = 6;  % horas

dydt = @(t,y) [8.63383E+08; %Crecimiento poblacional [mL^-1]
               (r.*I1./y(1))./(C + I1./y(1)) + b - d.*y(2); %Crecimiento protéico
               1.29998E+09; 
               (r.*I2./y(3))./(C + I2./y(3)) + b - d.*y(4); 
               4.87913E+07.*t*2 + 7.83194E+08;  
               (r.*I3./y(5))./(C + I3./y(5)) + b - d.*y(6); 
               1.62679E+08.*t*2 + 4.24091E+08; 
               (r.*I4./y(7))./(C + I4./y(7)) + b - d.*y(8)
               ];

[t, y] = ode45(dydt, 0:h:tlim, [p1;prot0;p2;prot0;p3;prot0;p4;prot0]);


linewidth = 2;

for con=1:8
    figure(con)%1)%
    plot(t, y(:,con),'LineWidth',linewidth)
    title('ConcentraciÃ³n')
    legend('Location','southeast')
    lgd = legend('show');
    lgd.FontSize = 20;
    xlabel('Time (a.u.)')
    ylabel('Concentration (a.u.)')
    label.FontSize = 20;
    grid on
    %hold on%
end
%figure(2)

