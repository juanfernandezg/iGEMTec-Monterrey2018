%% Species codification
% 1 --> I    (Cas1  protein)
% 2 --> II   (Cas2 protein)
% 3 --> X    (Cas proteins complex, 4 Cas1 + 2 Cas2)
% 4 --> RT   (Retrotranscriptase)
% 5 --> MSR
% 6 --> Y    (RT annexed to MSR while retrotranscribing)
% 7 --> MSD
% 8 --> A    (Storage machinery, X annexed to MSD)

%% Aditionally
% k_N  --> N specie reaction constant
% g_N  --> N specie degradation
% V_N  --> N specie basal production velocity
% C_N  --> N specie maximum concentration
% k__N --> N specie reverse reaction constant
% P    --> Promotor concentration (specifically IPTG)

%% Preguntar a Victor sobre modelado de la candidad de repeticiones que son
% posibles (Church, 3 veces). Muerte del sistema. Modelado de
% probabilidades

% tol = odeset('RelTol',1e-2); %Reduced tolerance (more precision)
h = 1;       %step
tlim = 180; %time limit
ilim = 15;   %insertion limit

[t, y] = ode45(@DESystem_v4, 0:h:tlim, [0;0;0;0;0;0;0;0]);

%% Poisson's lambda
lambda=cumsum(h*y(:,8));

%% Poisson mesh
k = 0:1:ilim;

%% pre-alloc
Ins = zeros(length(y(:,8)),length(k));

for i=1:length(k)
    Ins(:,i) = exp(-lambda(:,1)).*(lambda(:,1)).^k(i)./factorial(k(i));
end

figure(1)
mesh(Ins)
xlabel('Insertions')
ylabel('Time (a.u.)')
zlabel('Probability')
axes

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
% xlabel('Time (min)')
% ylabel('MSD concentration (arbitrary)')
% grid on

%% A
figure(8)
plot(t, y(:,8),'LineWidth',linewidth)
title('A vs Time')
legend('A','Location','southeast')
lgd = legend('show');
lgd.FontSize = 16;
xlabel('Time (a.u.)')
ylabel('A concentration (a.u.)')
grid on
