function [res] = aprox(x)
%% Quasi-steady state approximation simulation script
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

%% Parameters
% r,C,d,b Parameters
load('kRT.mat','kRT')
load('kWT.mat','kWT')

% Cas1 production
r_1  = 4/6.*kWT(1);
C_1  = kWT(2);
d_1  = kWT(3);
b_1  = 4/6.*kWT(4);

% Cas2 production
r_2  = 2/6.*kWT(1);
C_2  = kWT(2);
d_2  = kWT(3);
b_2  = 2/6.*kWT(4);

% RT production
r_RT  = kRT(1);
C_RT  = kRT(2);
d_RT  = kRT(3);
b_RT  = kRT(4);

% Message Processing (msr-msd + RT) reaction advancement constant 
k_M  = 1;

%% Differential Equations System solver
% tol = odeset('RelTol',1e-2); %Reduced tolerance (more precision)
format long
h    = .1;     %step size
tlim = 24;     %time limit
ilim =  5;     %insertion limit
I0   =   0.1;     %initial inductor concentration

a0 = x(1);
a1 = x(2);
a2 = x(3);
a3 = x(4);
a4 = x(5);
a5 = x(6);
a6 = x(7);
a7 = x(8);
a8 = x(9);
a9 = x(10);
a10 = x(11);
a11 = x(12);

%% List of Differential Equations    
dydt =  @(t,y) [0    ;%1D[Cas1]
    
        0    ;%2D[Cas2]
        
        0    ;%3D[RT]                                
                                                    
        0    ;%4D[UM]
        
        r_RT.*y(9)./(y(9)+C_RT) + b_RT   - k_M.*y(5)            ;%5D[MP]
        
        0    ;% 6 D[M]
        
        0    ;% 7 D[X]
        
        0    ;% 8 D[SM]
        
        % 9 D[I]                              
        %-2*(t-3).*exp(-(t-3).^2)       ;%Gaussian at t=3
        %2*cos(t).*sin(t)                        ;%Sin^2(t)  
        %0                                 ;%Constant IPTG
        a1 + 2*a2.*t.^1 + 3*a3.*t.^2 + 4*a4.*t.^3 + 5*a5.*t.^4;% + 6*a6.*t.^5 + 7*a7.*t.^6 + 8*a8.*t.^7 + 9*a9.*t.^8 + 10*a10.*t.^9 + 11*a11.*t.^10;%Taylor
        
        k_M.*y(5)   % 10 D[Ins]
        ];

[t, y] = ode45(dydt, 0:h:tlim, [0;0;0;0;0;0;0;0;a0;0]);  %Matlab solver

res = y(:,10);

end