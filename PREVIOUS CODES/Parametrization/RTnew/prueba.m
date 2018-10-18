%Población inicial
p1 = 1.14173E+09; % 0mM
p2 = 9.26913E+08; % 0.1mM
p3 = 8.49315E+08; % 0.5mM
p4 = 9.60488E+08; % 1.0mM
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

format long
h = .001;       %step
tlim = 6;  % horas

dydt = @(t,y) [-4.79040E+07.*t*2 + 4.35344E+08; %Crecimiento poblacional [mL^-1]
               (r.*I1./y(1))./(C + I1./y(1)) + b - d.*y(2); %Crecimiento protéico
               6.25191E+08; 
               (r.*I2./y(3))./(C + I2./y(3)) + b - d.*y(4); 
               -5.66984E+07.*t*2 + 6.78307E+08;  
               (r.*I3./y(5))./(C + I3./y(5)) + b - d.*y(6); 
               -1.51345E+07*3.*(t.^2) + 4.11483E+07.*t*2 + 6.07645E+08; 
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
end

