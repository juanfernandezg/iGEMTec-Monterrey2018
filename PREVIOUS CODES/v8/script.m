%% Versión 8, muy seguro de que está correcto
%% Species codification

% 1   --> 1      (Cas1 protein)
% 2   --> 2      (Cas2 protein)
% 3   --> RT     (RT protein)
% 4   --> UM     (Unprocessed message, msr-msd)
% 5   --> MP     (RT retrotranscribing msr-target sequence)
% 6   --> M      (Message, msDNA)
% 7   --> X      (Cas complex, 4 Cas1 + 2 Cas2)
% 8   --> SM     (Storage Machinery, X attached to msDNA)
% 9   --> I      (Inductor)

%% Constants nomenclature
% k_N  --> N specie reaction constant, i subindex for translation process
% d_N  --> N specie degradation constant
% r_N  --> N production rate
% C_N  --> N hill constant
% k__N --> N specie reverse reaction constant
% I    --> IPTG concentration

%% Preguntar a Victor sobre modelado de la candidad de repeticiones que son
% posibles (Church, 3 veces). Muerte del sistema. Modelado de
% probabilidades

% tol = odeset('RelTol',1e-2); %Reduced tolerance (more precision)
h = .01;       %step
tlim = 20; %time limit
ilim = 50;   %insertion limit
IPTG = 1;

[t, y] = ode45(@DESystem, 0:h:tlim, [zeros(8,1);IPTG;0]);

%% Poisson mesh
k = 0:1:ilim;

%% Surf
[tt,kk]=ndgrid(t,k);
Ins = exp(-y(:,10)).*y(:,10).^k./factorial(k);
figure(11)
surf(tt,kk,Ins,'EdgeColor','none')
linewidth = 2;

%% Gráficas 
% figure(2)
% plot(t, y(:,1),'LineWidth',linewidth)
% title('CasI vs Time')
% legend('Cas I','Location','southeast')
% lgd = legend('show');
% lgd.FontSize = 16;
% xlabel('Time (min)')
% ylabel('Cas I concentration (arbitrary)')
% grid on
% 
% figure(3)
% plot(t, y(:,2),'LineWidth',linewidth)
% title('Cas II vs Time')
% legend('Cas II','Location','southeast')
% lgd = legend('show');
% lgd.FontSize = 16;
% xlabel('Time (min)')
% ylabel('Cas II concentration (arbitrary)')
% grid on
% 
% figure(4)
% plot(t, y(:,3),'LineWidth',linewidth)
% title('Cas Complex vs Time')
% legend('Cas Complex','Location','southeast')
% lgd = legend('show');
% lgd.FontSize = 16;
% xlabel('Time (min)')
% ylabel('Cas Complex concentration (arbitrary)')
% grid on
% 
% figure(5)
% plot(t, y(:,4),'LineWidth',linewidth)
% title('RT vs Time')
% legend('RT','Location','southeast')
% lgd = legend('show');
% lgd.FontSize = 16;
% xlabel('Time (min)')
% ylabel('RT concentration (arbitrary)')
% grid on
% 
% figure(6)
% plot(t, y(:,5),'LineWidth',linewidth)
% title('MSR vs Time')
% legend('MSR','Location','southeast')
% lgd = legend('show');
% lgd.FontSize = 16;
% xlabel('Time (min)')
% ylabel('MSR concentration (arbitrary)')
% grid on
% 
% figure(7)
% plot(t, y(:,7),'LineWidth',linewidth)
% title('MSD vs Time')
% legend('MSD','Location','southeast')
% lgd = legend('show');
% lgd.FontSize = 16;
% % xlabel('Time (min)')
% % ylabel('MSD concentration (arbitrary)')
% % grid on
% 
% %% A
% figure(8)
% plot(t, y(:,24),'LineWidth',linewidth)
% title('A vs Time')
% legend('A','Location','southeast')
% lgd = legend('show');
% lgd.FontSize = 16;
% xlabel('Time (min)')
% ylabel('A concentration (arbitrary)')
% grid on


for con=1:10
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
    %hold on%
end
