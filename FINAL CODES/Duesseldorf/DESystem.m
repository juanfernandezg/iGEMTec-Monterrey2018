function dydt = DESystem(t,y)
% Species codification
% A --> E.Coli
% B --> S. Cerevisiae
% C --> S. Elongatus
% P --> Transport Protein


% Aditionally
% k_ij  --> constant of production by j of resource used by i
% c_i   --> consume constant of i
% Max_i --> Max concentration of ik__N --> -N specie constant
% p     --> Transport protein constant

% Parameters
c_A = rand(1);
c_B = rand(1);
c_C = rand(1);

b_A = rand(1);
b_B = rand(1);
b_C = rand(1);

A_max = rand(1);
B_max = rand(1);
C_max = rand(1);

k_AB = rand(1);
k_AC = rand(1);

k_BA = rand(1);
k_BC = rand(1);

k_CA = rand(1);
k_CB = rand(1);

% List of Differential Equations    
%{
dydt = [c_A*y(1)*y(2)*(max((1-y(1)/A_max),0) + b_A*y(1).^(1/2))                                ; % d[A]/dt
                 k_AB*y(3) + k_AC*y(5) - c_A*y(1)*(max((1-y(1)/A_max),0) + b_A*y(1).^(1/2))    ; % d[iA]/dt
        c_B*y(3)*y(4)*(max((1-y(2)/B_max),0) + b_B*y(2).^(1/2))                                ; % d[B]/dt
                 k_BA*y(1) + k_BC*y(5) - c_B*y(3)*(max((1-y(2)/B_max),0) + b_B*y(2).^(1/2))    ; % d[iB]/dt
        c_C*y(5)*y(6)*(max((1-y(3)/C_max),0) + b_C*y(3).^(1/2))                                ; % d[C]/dt
                 k_CA*y(1) + k_CB*y(3) - c_C*y(5)*(max((1-y(3)/C_max),0) + b_C*y(3).^(1/2))    ; % d[iC]/dt
];
%}
%{
dydt = [c_A*y(1)*y(2)*(max((1-y(1)/A_max),0) + b_A*y(1).^(2/3))                   ; % d[A]/dt
                 k_AB*y(3) + k_AC*y(5) - c_A*y(1)*y(2)*(max((1-y(1)/A_max),0) + b_A*y(1).^(2/3))      ; % d[iA]/dt
        c_B*y(3)*y(4)*(max((1-y(3)/B_max),0) + b_B*y(3).^(2/3))                   ; % d[B]/dt
                 k_BA*y(1) + k_BC*y(5) - c_B*y(3)*y(4)*(max((1-y(3)/B_max),0) + b_B*y(3).^(2/3))      ; % d[iB]/dt
        c_C*y(5)*y(6)*(max((1-y(5)/C_max),0) + b_C*y(5).^(2/3))                   ; % d[C]/dt
                 k_CA*y(1) + k_CB*y(3) - c_C*y(5)*y(6)*(max((1-y(5)/C_max),0) + b_C*y(5).^(2/3))     ; % d[iC]/dt
];
%}

dydt = [c_A*y(1)*y(2)*(1-y(1)/A_max)                                       ; % d[A]/dt
                 k_AB*y(3) + k_AC*y(5) - c_A*y(1)*y(2)*(1-y(1)/A_max)      ; % d[iA]/dt
        c_B*y(3)*y(4)*(1-y(3)/B_max)                                       ; % d[B]/dt
                 k_BA*y(1) + k_BC*y(5) - c_B*y(3)*y(4)*(1-y(3)/B_max)      ; % d[iB]/dt
        c_C*y(5)*y(6)*(1-y(5)/C_max)                                       ; % d[C]/dt
                 k_CA*y(1) + k_CB*y(3) - c_C*y(5)*y(6)*(1-y(5)/C_max)      ; % d[iC]/dt
];


end

