function dydt = DESystem(t,y)
% Species codification
% A --> E. Coli
% B --> S. Cerevisiae
% C --> S. Elongatus
% P --> Transport Protein (For a more descriptive model that considers a
% protein is responsible of nutrient transport from the medium into the
% bacteria

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

% System of Differential Equations
dydt = [c_A*y(1)*y(2)*(1-y(1)/A_max)                                       ; % d[A]/dt
                 k_AB*y(3) + k_AC*y(5) - c_A*y(1)*y(2)*(1-y(1)/A_max)      ; % d[iA]/dt
        c_B*y(3)*y(4)*(1-y(3)/B_max)                                       ; % d[B]/dt
                 k_BA*y(1) + k_BC*y(5) - c_B*y(3)*y(4)*(1-y(3)/B_max)      ; % d[iB]/dt
        c_C*y(5)*y(6)*(1-y(5)/C_max)                                       ; % d[C]/dt
                 k_CA*y(1) + k_CB*y(3) - c_C*y(5)*y(6)*(1-y(5)/C_max)        % d[iC]/dt
];
end

