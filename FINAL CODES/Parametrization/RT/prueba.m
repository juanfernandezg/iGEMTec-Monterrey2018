%Population (in case it was relevant for the model)
p1 = 1483727180; % 0mM
p2 = 1550066298; % 0.1mM
p3 = 1415314966; % 0.5mM
p4 = 1591528246; % 1.0mM

prot0 = [0.475165563;
0.723509934;
0.508278146;
0.706953642];
prot1=prot0(1); 
prot2=prot0(2); 
prot3=prot0(3); 
prot4=prot0(4); 

r=x(1);
C=x(2);
d=x(3);
b=x(4);
n=1;
di=0;

%Initial conditions
I1 = 0;
I2 = 0.1;
I3 = 0.5;
I4 = 1;

format long
h = .05;   %step size
tlim = 6;  %hours

dydt = @(t,y) [-4.79040E+07.*t*2 + 4.35344E+08; %Population growth [mL^-1]
               ((r.*I1.^n)./(C.^n + I1.^n) + b).*(1) - d.*y(2) - di.*y(1).*y(2); %Protein growth
               6.25191E+08; 
               ((r.*I2.^n)./(C.^n + I2.^n) + b).*(1) - d.*y(4) - di.*y(3).*y(4); 
               -5.66984E+07.*t*2 + 6.78307E+08;  
               ((r.*I3.^n)./(C.^n + I3.^n) + b).*(1) - d.*y(6) - di.*y(5).*y(6); 
               -1.51345E+07*3.*(t.^2) + 4.11483E+07.*t*2 + 6.07645E+08; 
               ((r.*I4.^n)./(C.^n + I4.^n) + b).*(1) - d.*y(8) - di.*y(7).*y(8)
               ];

[t, y] = ode45(dydt, 1:h:tlim, [p1;prot1;p2;prot2;p3;prot3;p4;prot4]);


linewidth = 2;

for con=1:8
    figure(con)%1)%
    plot(t, y(:,con),'LineWidth',linewidth)
    title('Concentración')
    legend('Location','southeast')
    lgd = legend('show');
    lgd.FontSize = 20;
    xlabel('Time (a.u.)')
    ylabel('Concentration (a.u.)')
    label.FontSize = 20;
    grid on
end

