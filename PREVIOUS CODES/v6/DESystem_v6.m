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

%% Constants nomenclature
% k_N  --> N specie reaction constant, i subindex for translation process
% d_N  --> N specie degradation constant
% r_N  --> N production rate
% C_N  --> N hill constant
% k__N --> N specie reverse reaction constant
% I    --> IPTG concentration
 
%% Parameters
xd = .1;

bRT = xd;
b1  = xd;
b2  = xd;

% IPTG concentration
I = xd;

% Cas1 production
r_1  = xd;
C_1  = xd;
d_1  = 0;

% Cas2 production
r_2  = xd;
C_2  = xd;
d_2  = 0;

% RT production
r_RT  = xd;
C_RT  = xd;
d_RT  = 0;

% Unprocessed Message (msr-msd) production
r_UM = xd;
d_UM = 0;

% Message Processing (msr-msd + RT)
k_MP  = xd;
d_MP  = 0;

% Message (msDNA) production rate
k_M  = xd;
d_M  = 0;

% Cas Complex 'X' (4 Cas1 + 2 Cas2)
k_X  = xd;
k__X = xd;
d_X = 0;

% Storage Machinery (X + msDNA) formation and insertion rates
k_SM  = xd;
k__SM = xd;
d_SM  = 0;
P_SM  = 1;

%% List of Differential Equations    
dydt = [4*b1 + r_1.*y(9)./(y(9)+C_1) - d_1.*y(1) - 4*k_X.*(y(1).^4).*(y(2).^2) + 4*k__X.*y(7)    ;%1D[Cas1]
        2*b2 + r_2.*y(9)./(y(9)+C_2) - d_2.*y(2) - 2*k_X.*(y(1).^4).*(y(2).^2) + 2*k__X.*y(7)    ;%2D[Cas2]
        
        bRT + r_RT.*y(9)./(y(9)+C_RT) - d_RT.*y(3) - k_MP.*y(3).*y(4) + d_MP.*y(5) + k_M.*y(5)   ;%3D[RT]
                                                    %Revisar si msr-msd se
                                                    %consume con la RT
        r_UM - d_UM.*y(4) - k_MP.*y(3).*y(4)                                                     ;% 4 D[UM], Ana informó que la producción del msr-ts-msd es constante sin inducción
        
        k_MP.*y(3).*y(4) - d_MP.*y(5) - k_M.*y(5)                                                ;% 5 D[MP]
        
        k_M.*y(5) - k_SM.*y(6).*y(7) + k__SM.*y(8) - d_M.*y(6)                                   ;% 6 D[M]
        
        k_X.*(y(1).^4).*(y(2).^2) - k__X.*y(7) - d_X.*y(7) - k_SM.*y(6).*y(7) + k__SM.*y(8) + P_SM.*y(8) ;% 7 D[X]
        
        k_SM.*y(6).*y(7) - k__SM.*y(8) - P_SM.*y(8) - k_SM.*y(8)                                         ;% 8 D[SM]
        
        % 9 D[Inductor]                              
        %-2*(t-10).*exp(-(t-10).^2)                  ;%Gaussiana en t=10
        %cos(t-10)                                      ;%Sin en t = 10
        0;
        ];
end