%% Complete Differential Equations System
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

%% Arbitrary Parameters, using values greater than 1 may result in long simulation times
ap = 1.5;

% Unprocessed Message (msr-msd) production
r_UM = ap;
C_UM = 0.01*ap;
d_UM = 0;
b_UM = ap;

% Message Processing (msr-msd + RT)
k_MP  = ap;
d_MP  = 0;

% Message (msDNA) production rate
k_M  = ap;
d_M  = 0;

% Cas Complex 'X' (4 Cas1 + 2 Cas2)
k_X  = ap;
k__X = 0;

% Storage Machinery (X + msDNA) formation and insertion rates
k_SM  = 0.9*ap;
d_SM  = 0;
P_SM  = 0.8*ap;

%% List of Differential Equations    
dydt = [b_1 + r_1.*y(9)./(y(9)+C_1) - d_1.*y(1) - 4*k_X.*(y(1).^4).*(y(2).^2) + 4*k__X.*y(7)     ;%1D[Cas1]
    
        b_2 + r_2.*y(9)./(y(9)+C_2) - d_2.*y(2) - 2*k_X.*(y(1).^4).*(y(2).^2) + 2*k__X.*y(7)     ;%2D[Cas2]        
        
        b_RT + r_RT.*y(9)./(y(9)+C_RT) - d_RT.*y(3) - k_MP.*y(3).*y(4) + d_MP.*y(5) + k_M.*y(5)  ;%3D[RT]                                                   
        
        b_UM + r_UM.*y(9)./(y(9)+C_UM) - d_UM.*y(4) - k_MP.*y(3).*y(4)                           ;%4D[UM]
        
        k_MP.*y(3).*y(4) - d_MP.*y(5) - k_M.*y(5)                                                ;%5D[MP]
        
        k_M.*y(5) - k_SM.*y(6).*y(7) - d_M.*y(6)                                                 ;%6D[M]
        
        k_X.*(y(1).^4).*(y(2).^2) - k__X.*y(7) - k_SM.*y(6).*y(7) + d_SM.*y(8) + P_SM.*y(8)      ;%7D[X]
        
        k_SM.*y(6).*y(7) - d_SM.*y(8) - P_SM.*y(8)                                               ;%8D[SM]
        
        % 9 D[Inductor]                              
        %-2*(t-10).*exp(-(t-10).^2)                  ;%Gaussiana en t=10
        %cos(t-10)                                   ;%Sin en t = 10
        0;
        
        P_SM.*y(8)  %10D[Ins]
];
end