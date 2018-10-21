%% Quasi-steady state approximation Differential Equations System
function dydt = DESystem(t,y)
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

%% List of Differential Equations    
dydt = [0    ;%1D[Cas1]
    
        0    ;%2D[Cas2]
        
        0    ;%3D[RT]                                
                                                    
        0    ;%4D[UM]
        
        r_RT.*y(9)./(y(9)+C_RT) + b_RT   - k_M.*y(5)            ;%5D[MP]
        
        0    ;% 6 D[M]
        
        0    ;% 7 D[X]
        
        0    ;% 8 D[SM]
        
        % 9 D[I]                              
        %-2*(t-3).*exp(-(t-3).^2)       ;%Gaussian at t=3
        %2*cos(t/8).*sin(t/8)                        ;%Sin^2(t)  
        %0                                 ;%Constant IPTG
        2*t - (t.^2)/10;%
        
        k_M.*y(5)   % 10 D[Ins]
        ];
end