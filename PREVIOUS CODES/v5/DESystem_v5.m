function dydt = DESystem(t,y)
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
% 10   --> m2    (mRNA for Cas2)
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
% d_N    --> N specie degradation constant
% dr_N   --> N mRNA ribosome drop-off constant
% r_N    --> N mRNA production rate
% C_N    --> N specie maximum concentration
% k__N   --> N specie reverse reaction constant
% IPTG   --> IPTG concentration
 
%% Parameters
xd = .1;

% IPTG concentration
I = xd;

% Cas1 production
r_1    = xd;
C_1    = xd;
k_10   = xd;
k_11   = xd;
k_12   = xd;
k_13   = xd;
k_14   = xd;
d_m01  = xd;
d_b01  = xd;
d_m1   = xd;
d_b1   = xd;
d_R1   = xd;
dr_1   = xd;
d_Cas1 = xd;

% Cas2 production
r_2    = xd;
C_2    = xd;
k_20   = xd;
k_21   = xd;
k_22   = xd;
k_23   = xd;
k_24   = xd;
d_m02  = xd;
d_b02  = xd;
d_m2   = xd;
d_b2   = xd;
d_R2   = xd;
dr_2   = xd;
d_Cas2 = xd;

% RT production
r_RT   = xd;
C_RT   = xd;
k_RT0  = xd;
k_RT1  = xd;
k_RT2  = xd;
k_RT3  = xd;
k_RT4  = xd;
d_m0RT = xd;
d_b0RT = xd;
d_mRT  = xd;
d_bRT  = xd;
d_RRT  = xd;
dr_RT  = xd;
d_RT   = xd;

% Unprocessed Message (msr-msd) production
r_UM = xd;
C_UM = xd;
d_UM = xd;

% Message Processing (msr-msd + RT) rate
k_MP  = xd;
k__MP = xd;

% Message (msDNA) production rate
k_M  = xd;
k__M = xd;

% Cas Complex 'X' (4 Cas1 + 2 Cas2)
k_X  = xd;
k__X = xd;

bRT=xd;

b1=xd;
b2=xd;

% Storage Machinery (X + msDNA) formation and insertion rates
k_SM  = xd;
k__SM = xd;
P_SM  = 1;


%% List of Differential Equations    
dydt = [-y(1)                                                              ;% D[IPTG]
    
   b1+  r_1.*y(1)./(y(1)+C_1)   - (d_m01 + k_10).*y(2)                          ;% D[m01]
        k_10.*y(2)              - (d_b01 + k_11).*y(3)                          ;% D[b01]
        k_11.*y(3) + k_13.*y(5) - (d_m1  + k_12).*y(4)                          ;% D[m1]
        k_12.*y(4)              - (d_b1  + k_13).*y(5)                          ;% D[b1]
        k_11.*y(3) + k_13.*y(5) - (d_R1 + dr_1 + k_14).*y(6)                    ;% D[R1]
        k_14.*y(6) + k__X.*y(23) -  d_Cas1.*y(7) - k_X.*(y(7).^4).*(y(13).^2)   ;% D[Cas1]
        
   b2+  r_2.*y(1)./(y(1)+C_2)    - (d_m02 + k_20).*y(8)                         ;% D[m02]
        k_20.*y(8)               - (d_b02 + k_21).*y(9)                         ;% D[b02]
        k_21.*y(9) + k_23.*y(11) - (d_m2  + k_22).*y(10)                        ;% D[m2]
        k_22.*y(10)              - (d_b2  + k_23).*y(11)                        ;% D[b2]
        k_21.*y(9) + k_23.*y(11) - (d_R2 + dr_2 + k_24).*y(12)                  ;% D[R2]
        k_24.*y(12) + k__X.*y(23) -  d_Cas2.*y(13) - k_X.*(y(7).^4).*(y(13).^2) ;% D[Cas2]
                                   
   bRT+ r_RT.*y(1)./(y(1)+C_RT)     - (d_m0RT + k_RT0).*y(14)                       ;% D[m0RT]
        k_RT0.*y(14)                - (d_b0RT + k_RT1).*y(15)                       ;% D[b0RT]
        k_RT1.*y(15) + k_RT3.*y(17) - (d_mRT  + k_RT2).*y(16)                       ;% D[mRT]
        k_RT2.*y(16)                - (d_bRT  + k_RT3).*y(17)                       ;% D[bRT]
        k_RT1.*y(15) + k_RT3.*y(17) - (d_RRT + dr_RT + k_RT4).*y(18)                ;% D[RRT]
        k_RT4.*y(18) + k_M.*y(21) + k__MP.*y(21) - d_RT.*y(19) - k_MP.*y(19).*y(20) ;% D[RT]
        
        r_UM - d_UM.*y(20) - k_MP.*y(19).*y(20) + k__MP.*y(21) ;% D[UM], Ana informó que la producción del msr-ts-msd es constante sin inducción
        
        k_MP.*y(19).*y(20) - k__MP.*y(21) - k_M.*y(21)         ;% D[MP]
        
        k_M.*y(21) - k_SM.*y(22).*y(23) + k__SM.*y(24)         ;% D[M]
        
        k_X.*(y(7).^4).*(y(13).^2) - k__X.*y(23) - k_SM.*y(22).*y(23) + P_SM.*y(24) ;% D[X]
        
        k_SM.*y(22).*y(23) - k__SM.*y(24) - P_SM.*y(24)        ;% D[SM]
];
end