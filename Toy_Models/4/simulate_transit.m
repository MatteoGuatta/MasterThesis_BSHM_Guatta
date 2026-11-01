function [acc_sensors, t_target] = simulate_transit(cfg, V, P_val, damage_info)

%cfg: File di configurazione contenente tutti i parametri, geometrici e
%operativi 
%V: Velocità del convoglio
%P_val: Intesità di carico


%% Introduzione dello scenario di danno 

EJ_curr = cfg.EJ ;
EA_curr = cfg.EA; 

if damage_info.is_damaged == true
    el = damage_info.element;
    EJ_curr(el) = EJ_curr(el) * (1 - damage_info.severity);
    EA_curr(el) = EA_curr(el) * (1 - damage_info.severity);
end


%% Assemblaggio delle matrici M,C,K

[M, K, C] = assem(cfg.incidenze, cfg.l, cfg.m, EA_curr, EJ_curr, ...
                  cfg.gamma, cfg.idb, cfg.conc_data, cfg.rayleigh_ab);


    
% Partizione sui singoli gradi di libertà
MFF = M(1:cfg.ndof, 1:cfg.ndof);
CFF = C(1:cfg.ndof, 1:cfg.ndof);
KFF = K(1:cfg.ndof, 1:cfg.ndof);


%% Simulazione

% Posizione iniziale dei carichi
pos_P_in = [];
x_curr = 0;

for iv = 1:cfg.n_vagoni 
    
    % Posizionamento dei primi due assi del vagone iv
    pos_P_in = [pos_P_in; -x_curr; -(x_curr + cfg.D1)];

    % Avanzamento al vagone successivo
    x_curr = x_curr + cfg.D1 + cfg.D2;

end

n_P = length(pos_P_in);
P_vec = P_val * ones(n_P,1);

% Intervallo temporale di integrazione
deck_el      = cfg.deck_elements(:);
x_start_deck = cfg.posiz(deck_el, 1);
l_deck       = cfg.l(deck_el); l_deck = l_deck(:);
gamma_deck   = cfg.gamma(deck_el); gamma_deck = gamma_deck(:);
x_end_deck   = x_start_deck + l_deck .* cos(gamma_deck);

x_deck_start = min([x_start_deck; x_end_deck]);
x_deck_end   = max([x_start_deck; x_end_deck]);

% Durata della simulazione: usiamo la durata fissa cfg.t_total (15 s)
% per garantire che tutte le simulazioni abbiano esattamente la stessa dimensione.
t_exit = (x_deck_end - min(pos_P_in)) / V; % Istante teorico di uscita dell'ultimo asse

dt = cfg.dt_newmark;                                 
n_steps = ceil(cfg.t_total / dt);
t = linspace(0, cfg.t_total, n_steps + 1)';

% Funzione forzante dinamica
F_fun = @(tt) Forze_carichi_viaggianti( ...
    tt, cfg.ndof, cfg.xy, cfg.incidenze, cfg.l, cfg.gamma, cfg.posiz, cfg.deck_elements, pos_P_in, V, P_vec);

%% Integrazione diretta con algoritmo di Newmark
beta_N  = 1/4; % Variante accelerazione media costante (incondizionatamente stabile)
gamma_N = 1/2;

q0     = zeros(cfg.ndof, 1);
q_dot0 = zeros(cfg.ndof, 1);

[q, q_dot, q_ddot] = Newmark_accelerazioni( ...
    MFF, CFF, KFF, F_fun, t, q0, q_dot0, beta_N, gamma_N, false);

%% Estrazione dei dati accelerometrici relativi ai sensori

% 1. Preallocazione 
X = zeros(cfg.n_sensors, size(q_ddot, 2));

% 2. Ciclo sui 9 sensori
for i = 1:cfg.n_sensors
    nodo = cfg.sensor_nodes(i);
    dof_w = cfg.idb(nodo, 2);
    X(i, :) = q_ddot(dof_w, :);
end

%% Ricampionamento alla frequenza target dei sensori (200 Hz -> esattamente 3000 punti)
t_target = (0:cfg.pts_total-1)' / cfg.fs_sensor;    % Vettore tempo da 0 a 14.995 s
acc_sensors = interp1(t, X', t_target, 'linear')'; % [9 x 3000]

%% Aggiunta rumore di misura del sensore (5% AWGN)
if isfield(cfg, 'noise_sensor') && cfg.noise_sensor > 0
    for is = 1:cfg.n_sensors
        sig_std = std(acc_sensors(is, :));
        noise = cfg.noise_sensor * sig_std * randn(1, cfg.pts_total);
        acc_sensors(is, :) = acc_sensors(is, :) + noise;
    end
end

end







    