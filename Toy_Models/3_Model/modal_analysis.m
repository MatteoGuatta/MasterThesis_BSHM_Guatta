function [freq, Phi, zeta, modal_par_list_table] = modal_analysis(M, K, C, moi)
% MODAL_ANALYSIS Esegue l'analisi modale del sistema dinamico ridotto sui GDL liberi
%
% Input:
%   M, K, C : matrici strutturali [ndof x ndof] (sui soli GDL liberi)
%   moi     : numero di modi di vibrare da estrarre
%
% Output:
%   freq                 : frequenze proprie [Hz]
%   Phi                  : forme modali normalizzate a 1
%   zeta                 : smorzamenti adimensionati dei modi
%   modal_par_list_table : tabella riassuntiva dei parametri modali

Ndof_f = size(M, 1);
if moi > Ndof_f
    error(['Numero di modi richiesti (', num2str(moi), ') maggiore dei GDL disponibili (', num2str(Ndof_f), ').']);
end

modal_par_list = zeros(moi, 7);

[Phi, Lambda] = eig(K, M);

omega2 = diag(Lambda);
omega  = sqrt(abs(omega2));

[omega, idx] = sort(omega);
Phi = Phi(:, idx);

% Normalizzazione delle forme modali al massimo valore assoluto unitario
for ik = 1:moi
    max_val = max(abs(Phi(:, ik)));
    if max_val > 0
        Phi(:, ik) = Phi(:, ik) ./ max_val;
    end
end

C_mod = Phi' * C * Phi;
M_mod = Phi' * M * Phi;
K_mod = Phi' * K * Phi;

Phi  = Phi(:, 1:moi);
freq = omega(1:moi) / (2 * pi);

c_modal_tot = diag(C_mod);
m_modal_tot = diag(M_mod);
k_modal_tot = diag(K_mod);

zeta = zeros(moi, 1);
for ik = 1:moi
    if m_modal_tot(ik) > 0 && omega(ik) > 0
        zeta(ik) = c_modal_tot(ik) / (2 * m_modal_tot(ik) * omega(ik));
    end
end

modal_par_list(:, 1) = (1:moi)';
modal_par_list(:, 2) = m_modal_tot(1:moi);
modal_par_list(:, 3) = k_modal_tot(1:moi);
modal_par_list(:, 4) = c_modal_tot(1:moi);
modal_par_list(:, 5) = omega(1:moi);
modal_par_list(:, 6) = freq(1:moi);
modal_par_list(:, 7) = zeta(1:moi);

modal_par_list_table = array2table(modal_par_list, ...
    'VariableNames', {'Mode', 'm_modal_kg', 'k_modal_Nm', 'c_modal_Nsm', 'omega_rad_s', 'freq_Hz', 'zeta'});

end
