function dydt = DESystem_v4(t,y)
%% Species codification
% 1 --> I    (Cas1  protein)
% 2 --> II   (Cas2 protein)
% 3 --> X    (Cas proteins complex, 4 Cas1 + 2 Cas2)
% 4 --> RT   (Retrotranscriptase
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
 
%% Parameters
P    = 0.1;

V_I  = 0.1;
C_I  = 0.1;
g_I  = 0.1;

V_II = V_I;
C_II = 0.1;
g_II = 0.1;

k_X  = 0.1;
k__X = 0.1;

V_RT = 0.1;
C_RT = 0.1;
g_RT = 0.1;

V_MSR  = 0.1;
C_MSR  = 0.1;
g_MSR  = 0.1;

k_Y   = 0.1;
k__Y  = 0.1;

k_MSD  = 0.1;
g_MSD  = 0.1;

k_A  = 0.1;
k__A = 0.1;

%% List of Differential Equations    
dydt = [r_1.*y(9)./(y(9)+C_1) - d_1.*y(1) - 4*k_X.*(y(1).^4).*(y(2).^2) + 4*k__X.*y(7)    ;%D[Cas1]
        r_2.*y(9)./(y(9)+C_2) - d_2.*y(2) - 2*k_X.*(y(1).^4).*(y(2).^2) + 2*k__X.*y(7)    ;%D[Cas2]
        
        r_RT.*y(9)./(y(9)+C_RT) - d_RT.*y(3) - k_MP.*y(3).*y(4) + d_MP.*y(5) + k_M.*y(5)  ;%D[RT]
        
        r_UM - d_UM.*y(4) - k_MP.*y(3).*y(4) ;% D[UM], Ana informó que la producción del msr-ts-msd es constante sin inducción
        
        k_MP.*y(3).*y(4) - d_MP.*y(5) - k_M.*y(5)         ;% D[MP]
        
        k_M.*y(5) - k_SM.*y(6).*y(7) + k__SM.*y(8)         ;% D[M]
        
        k_X.*(y(1).^4).*(y(2).^2) - k__X.*y(7) - k_SM.*y(6).*y(7) + k__SM.*y(8) + P_SM.*y(8) ;% D[X]
        
        k_SM.*y(6).*y(7) - k__SM.*y(8) - P_SM.*y(8)                                          ;% D[SM]
        
        0                   ;%Inductor
        ];
end