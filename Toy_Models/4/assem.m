function [M, K, C] = assem(incidenze, l, m, EA, EJ, gamma, idb, conc_data, rayleigh_ab)
% ASSEM Assembla le matrici globali di massa, rigidezza e smorzamento
% Include i contributi delle travi 2D (assiale + flessionale) e degli
% elementi concentrati (molle al suolo, masse/inerzie, smorzatori).
%
% Input:
%   incidenze   - [n_el x 6] Matrice delle incidenze dei gradi di libertà
%   l           - [1 x n_el] Lunghezze degli elementi
%   m           - [1 x n_el] Masse lineari [kg/m]
%   EA          - [1 x n_el] Rigidezze assiali [N]
%   EJ          - [1 x n_el] Rigidezze flessionali [N*m^2]
%   gamma       - [1 x n_el] Angoli di orientazione rispetto all'asse X globale [rad]
%   idb         - [nnod x 3] Matrice dei gradi di libertà
%   conc_data   - (Opzionale) Struct contenente elementi concentrati:
%                   .springs: [nodo, kx, ky, ktheta]
%                   .masses:  [nodo, mx, my, Jtheta]
%                   .dampers: [nodo, cx, cy, ctheta]
%   rayleigh_ab - (Opzionale) Vettore [alpha, beta] per lo smorzamento di Rayleigh
%
% Output:
%   M           - Matrice di massa globale [n_dof_tot x n_dof_tot]
%   K           - Matrice di rigidezza globale [n_dof_tot x n_dof_tot]
%   C           - Matrice di smorzamento globale [n_dof_tot x n_dof_tot]

%% 1. Controlli di consistenza
n_el = size(incidenze, 1);
if size(incidenze, 2) ~= 6
    error('Errore: la matrice incidenze deve avere 6 colonne (3 GDL per nodo).');
end
if length(l) ~= n_el || length(m) ~= n_el || length(EA) ~= n_el || length(EJ) ~= n_el || length(gamma) ~= n_el
    error('Errore: dimensioni dei vettori delle proprietà non coerenti con il numero di elementi.');
end

n_dof_tot = max(idb(:));

%% 2. Inizializzazione matrici
M_beam = zeros(n_dof_tot, n_dof_tot);
K_beam = zeros(n_dof_tot, n_dof_tot);
C_conc = zeros(n_dof_tot, n_dof_tot);

%% 3. Assemblaggio contributi travi (distribuiti)
for k = 1:n_el
    [mG, kG] = el_tra(l(k), m(k), EA(k), EJ(k), gamma(k));
    for iri = 1:6
        for ico = 1:6
            i1 = incidenze(k, iri);
            i2 = incidenze(k, ico);
            M_beam(i1, i2) = M_beam(i1, i2) + mG(iri, ico);
            K_beam(i1, i2) = K_beam(i1, i2) + kG(iri, ico);
        end
    end
end

M = M_beam;
K = K_beam;

%% 4. Contributi di elementi concentrati (se forniti)
if nargin >= 8 && ~isempty(conc_data)
    % A. Molle concentrate a terra
    if isfield(conc_data, 'springs') && ~isempty(conc_data.springs)
        for s = 1:size(conc_data.springs, 1)
            node = conc_data.springs(s, 1);
            k_val = conc_data.springs(s, 2:4);
            for d = 1:3
                i_dof = idb(node, d);
                if i_dof > 0
                    K(i_dof, i_dof) = K(i_dof, i_dof) + k_val(d);
                end
            end
        end
    end

    % B. Masse e inerzie concentrate nodali
    if isfield(conc_data, 'masses') && ~isempty(conc_data.masses)
        for s = 1:size(conc_data.masses, 1)
            node = conc_data.masses(s, 1);
            m_val = conc_data.masses(s, 2:4);
            for d = 1:3
                i_dof = idb(node, d);
                if i_dof > 0
                    M(i_dof, i_dof) = M(i_dof, i_dof) + m_val(d);
                end
            end
        end
    end

    % C. Smorzatori concentrati
    if isfield(conc_data, 'dampers') && ~isempty(conc_data.dampers)
        for s = 1:size(conc_data.dampers, 1)
            node = conc_data.dampers(s, 1);
            c_val = conc_data.dampers(s, 2:4);
            for d = 1:3
                i_dof = idb(node, d);
                if i_dof > 0
                    C_conc(i_dof, i_dof) = C_conc(i_dof, i_dof) + c_val(d);
                end
            end
        end
    end
end

%% 5. Assemblaggio matrice di smorzamento C
if nargin >= 9 && ~isempty(rayleigh_ab)
    alpha = rayleigh_ab(1);
    beta  = rayleigh_ab(2);
    C = alpha * M_beam + beta * K_beam + C_conc;
else
    C = C_conc;
end

end
