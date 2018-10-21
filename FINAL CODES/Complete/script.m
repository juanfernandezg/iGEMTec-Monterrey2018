%% Complete Differential Equations System script
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
% 10  --> Ins    (Insertion)

%% Constants nomenclature
% k_N  --> N species reaction constant
% d_N  --> N species degradation constant
% r_N  --> N species production rate
% C_N  --> N species Hill equation constant
% k__N --> N species reverse reaction constant
% I0   --> Initial inductor concentration

%% Differential Equations System solver
% tol = odeset('RelTol',1e-2); %Reduced tolerance (more precision)
format long
h    = 1;     %step size
tlim = 400;     %time limit
ilim =  15;     %insertion limit
I0   =   0.1;     %initial inductor concentration

[t, y] = ode45(@DESystem, 0:h:tlim, [0;0;0;0;0;0;0;0;I0;0]);  %Matlab solver

%% Insertion probabilities distribution
k = 0:1:ilim;            %Vector space for k insertions
[tt,kk]=ndgrid(t,k);     %Vector space for k insertions at t time
Ins = exp(-y(:,10)).*y(:,10).^k./factorial(k);    %Poisson distribution
figure(11)
surf(tt,kk,Ins,'EdgeColor','none')     %Plot surface

%% Graphs
linewidth = 2;
for con=1:10
    figure(con)
    plot(t, y(:,con),'LineWidth',linewidth)
    title('Concentration')
    legend('Location','southeast')
    lgd = legend('show');
    lgd.FontSize = 20;
    xlabel('Time (hours)')
    ylabel('Concentration (mM)')
    label.FontSize = 20;
    grid on
    %hold on   %To graph them together
end
figure(11)
