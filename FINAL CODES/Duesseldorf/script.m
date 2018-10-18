% Species codification
% A --> E. Coli
% B --> S. Cerevisiae
% C --> S. Elongatus

%Reduced tolerance
%tol = odeset('RelTol',1e-2);

[t, y] = ode45(@DESystem, 0:.1:10, [rand(1);rand(1);rand(1);rand(1);rand(1);rand(1)]);

linewidth = 2;

figure(1)
plot(t, y(:,1),'LineWidth',linewidth)
title('E.Coli vs Time')
legend('E.Coli','Location','southeast')
lgd = legend('show');
lgd.FontSize = 16;
xlabel('Time (arbitrary)')
ylabel('E.Coli population (arbitrary)')
grid on

figure(2)
plot(t, y(:,3),'LineWidth',linewidth)
title('S. Cerevisiae vs Time')
legend('S. Cerevisiae','Location','southeast')
lgd = legend('show');
lgd.FontSize = 16;
xlabel('Time (arbitrary)')
ylabel('S. Cerevisiae population (arbitrary)')
grid on

figure(3)
plot(t, y(:,5),'LineWidth',linewidth)
title('S. Elongatus vs Time')
legend('S. Elongatus','Location','southeast')
lgd = legend('show');
lgd.FontSize = 16;
xlabel('Time (arbitrary)')
ylabel('S. Elongatus population (arbitrary)')
grid on
