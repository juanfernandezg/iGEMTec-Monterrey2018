%% Species codification

% 1   --> I      (IPTG concentration)

% 2   --> m01    (mRNA newborn for Cas1)
% 3   --> b01    (Bound mRNA newborn for Cas1)
% 4   --> m1     (mRNA for Cas1)
% 5   --> b1     (Bound mRNS for Cas1)
% 6   --> R1     (Active ribosome for Cas1)
% 7   --> Cas1   (Cas1 protein)

% 8   --> m02    (mRNA newborn for Cas2)
% 9   --> b02    (Bound mRNA newborn for Cas2)
% 10  --> m2    (mRNA for Cas2)
% 11  --> b2     (Bound mRNS for Cas2)
% 12  --> R2     (Active ribosome for Cas2)
% 13  --> Cas2   (Cas2 protein)

% 14  --> m0RT   (mRNA newborn for RT)
% 15  --> b0RT   (Bound mRNA newborn for RT)
% 16  --> mRT    (mRNA for RT)
% 17  --> bRT    (Bound mRNS for RT)
% 18  --> RRT    (Active ribosome for RT)
% 19  --> RT     (RT protein)

% 20  --> UM     (Unprocessed message, msr-msd)
% 21  --> MP     (RT retrotranscribing msr-target sequence)
% 22  --> M      (Message, msDNA)
% 23  --> X      (Cas complex, 4 Cas1 + 2 Cas2)
% 24  --> SM     (Storage Machinery, X attached to msDNA)

%% Constants nomenclature
% k_Ni  --> N specie reaction constant, i subindex for translation process
% d_N   --> N specie degradation constant
% dr_N  --> N mRNA ribosome drop-off constant
% r_N   --> N mRNA production rate
% C_N   --> N specie maximum concentration
% k__N  --> N specie reverse reaction constant
% IPTG  --> IPTG concentration
 

% Preguntar a Victor sobre modelado de la candidad de repeticiones que son
% posibles (Church, 3 veces). Muerte del sistema. Modelado de
% probabilidades



%% Preguntar a Victor sobre modelado de la candidad de repeticiones que son
% posibles (Church, 3 veces). Muerte del sistema. Modelado de
% probabilidades

% tol = odeset('RelTol',1e-2); %Reduced tolerance (more precision)
h = .1;       %step
tlim = 10000; %time limit
ilim = 100;   %insertion limit

[t, y] = ode45(@DESystem, 0:h:tlim, [10;zeros(23,1)]);

% Poisson's lambda
lambda=cumsum(h*y(:,24));

%% Poisson mesh
k = 0:1:ilim;

%% pre-alloc
Ins = zeros(length(y(:,24)),length(k));

for i=1:length(k)
    Ins(:,i) = exp(-lambda(:,1)).*(lambda(:,1)).^k(i)./factorial(k(i));
end

figure(25)
mesh(Ins)

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


for con=1:24
    figure(con)
    plot(t, y(:,con),'LineWidth',linewidth)
    title('Concentración')
    legend('Location','southeast')
    lgd = legend('show');
    lgd.FontSize = 16;
    xlabel('Time (min)')
    ylabel('A concentration (arbitrary)')
    grid on
end