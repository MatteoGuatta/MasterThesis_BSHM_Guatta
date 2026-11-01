function [mG, kG] = el_tra(l, m, EA, EJ, gamma)
% EL_TRA Matrici di massa e rigidezza dell'elemento trave 2D (6 GDL)
% Comprende deformazione assiale (EA) e flessione Euler-Bernoulli (EJ),
% ruotate nel sistema di riferimento globale mediante l'angolo gamma.
%
% GDL locali: [u1; w1; theta1; u2; w2; theta2]

% Matrice di massa locale consistente (assiale + flessionale)
mL = m * l * [ ...
    1/3,       0,          0,           1/6,      0,          0;
    0,         13/35,      11*l/210,    0,        9/70,      -13*l/420;
    0,         11*l/210,   l^2/105,     0,        13*l/420,  -l^2/140;
    1/6,       0,          0,           1/3,      0,          0;
    0,         9/70,       13*l/420,    0,        13/35,     -11*l/210;
    0,        -13*l/420,  -l^2/140,     0,       -11*l/210,   l^2/105  ];

% Matrice di rigidezza locale: contributo assiale
kL_ax = (EA / l) * [ ...
     1,  0,  0, -1,  0,  0;
     0,  0,  0,  0,  0,  0;
     0,  0,  0,  0,  0,  0;
    -1,  0,  0,  1,  0,  0;
     0,  0,  0,  0,  0,  0;
     0,  0,  0,  0,  0,  0 ];

% Matrice di rigidezza locale: contributo flessionale Euler-Bernoulli
kL_fl = EJ * [ ...
     0,      0,          0,       0,      0,          0;
     0,  12/l^3,     6/l^2,       0, -12/l^3,     6/l^2;
     0,   6/l^2,       4/l,       0,  -6/l^2,       2/l;
     0,      0,          0,       0,      0,          0;
     0, -12/l^3,    -6/l^2,       0,  12/l^3,    -6/l^2;
     0,   6/l^2,       2/l,       0,  -6/l^2,       4/l  ];

kL = kL_ax + kL_fl;

% Matrice di rotazione da riferimento locale a globale (3x3 per nodo)
lambda = [ cos(gamma),  sin(gamma), 0;
          -sin(gamma),  cos(gamma), 0;
           0,           0,          1 ];

% Matrice di rotazione per l'elemento a 6 GDL
Lambda = blkdiag(lambda, lambda);

% Trasformazione nel riferimento globale
mG = Lambda' * mL * Lambda;
kG = Lambda' * kL * Lambda;
end
