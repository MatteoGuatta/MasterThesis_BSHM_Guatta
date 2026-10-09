function [q, q_dot, q_ddot] = Newmark_accelerazioni( ...
    M, C, K, F_fun, t, q0, q_dot0, beta_N, gamma_N)
% NEWMARK_ACCELERAZIONI Risolve le equazioni del moto con algoritmo di Newmark implicito
% Restituisce storie temporali di spostamento, velocità e accelerazione.
%
% Input:
%   M, C, K : matrici strutturali [ndof x ndof]
%   F_fun   : funzione anonima @(tt) F(tt) che restituisce il vettore delle forze al tempo tt
%   t       : vettore degli istanti temporali di integrazione [nt x 1]
%   q0      : vettore spostamenti iniziali [ndof x 1]
%   q_dot0  : vettore velocità iniziali [ndof x 1]
%   beta_N  : parametro Newmark (es. 1/4 per accelerazione media costante)
%   gamma_N : parametro Newmark (es. 1/2)
%
% Output:
%   q       : [ndof x nt] spostamenti
%   q_dot   : [ndof x nt] velocità
%   q_ddot  : [ndof x nt] accelerazioni

dt = t(2) - t(1);
kd = size(M, 1);
nt = numel(t);

%% Inizializzazione vettori di stato
q      = zeros(kd, nt);
q_dot  = zeros(kd, nt);
q_ddot = zeros(kd, nt);

q(:, 1)     = q0(:);
q_dot(:, 1) = q_dot0(:);

% Accelerazione iniziale dall'equilibrio dinamico
q_ddot(:, 1) = M \ (F_fun(t(1)) - C * q_dot(:, 1) - K * q(:, 1));

%% Matrice dinamica efficace e fattorizzazione di Cholesky
A_eff = M + gamma_N * dt * C + beta_N * dt^2 * K;
L_eff = chol(A_eff, 'lower');

%% Integrazione temporale passo-passo
prossima_percentuale = 10;
fprintf('Simulation status: 0%%\n');

for n = 1:nt-1
    % Predittori
    q_pred = q(:, n) + dt * q_dot(:, n) + dt^2 * (0.5 - beta_N) * q_ddot(:, n);
    v_pred = q_dot(:, n) + dt * (1 - gamma_N) * q_ddot(:, n);

    % Forze al passo successivo t(n+1)
    F_next = F_fun(t(n + 1));

    % Risoluzione per l'accelerazione al passo successivo
    rhs = F_next - C * v_pred - K * q_pred;
    q_ddot(:, n + 1) = L_eff' \ (L_eff \ rhs);

    % Correttori (aggiornamento spostamento e velocità)
    q(:, n + 1)     = q_pred + beta_N * dt^2 * q_ddot(:, n + 1);
    q_dot(:, n + 1) = v_pred + gamma_N * dt * q_ddot(:, n + 1);

    % Avanzamento percentuale a terminale
    percentuale = floor(100 * n / (nt - 1));
    if percentuale >= prossima_percentuale
        fprintf('Simulation status: %d%%\n', percentuale);
        prossima_percentuale = prossima_percentuale + 10;
    end
end
fprintf('Simulation status: 100%% - Completed.\n');

end
