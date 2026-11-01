function [acc_sensors, t_target] = simulate_transit(cfg, struct_model, V, P_val, damage_info, apply_noise)
% SIMULATE_TRANSIT Esegue la simulazione dinamica FEM di un singolo transito
% ferroviario sul ponte, estraendo i segnali dei sensori virtuali.
%
% Input:
%   cfg          : struct generata da config_benchmark()
%   struct_model : struct con dati geometrici e topologici da loadstructure()
%   V            : velocità del convoglio [m/s]
%   P_val        : carico per asse [N]
%   damage_info  : struct con campi:
%                    .is_damaged (true/false)
%                    .element    (indice dell'elemento danneggiato, es. 15)
%                    .severity   (frazione di riduzione di EJ, es. 0.25 per -25%)
%   apply_noise  : (Opzionale, default=true) se true aggiunge rumore di misura e materiale
%
% Output:
%   acc_sensors  : [9 x 3000] matrice delle accelerazioni verticali sui 9 sensori [m/s^2]
%   t_target     : [3000 x 1] vettore temporale campionato a fs = 200 Hz [s]

    if nargin < 6 || isempty(apply_noise)
        apply_noise = true;
    end
    if nargin < 5 || isempty(damage_info)
        damage_info.is_damaged = false;
        damage_info.element    = [];
        damage_info.severity   = 0;
    end

    %% 1. Inizializzazione proprietà strutturali con eventuale incertezza naturale
    EJ_curr = struct_model.EJ_nom;
    EA_curr = struct_model.EA_nom;
    m_curr  = struct_model.m_nom;
    
    if apply_noise && isfield(cfg, 'sigma_E_mat') && cfg.sigma_E_mat > 0
        EJ_curr = EJ_curr .* (1 + cfg.sigma_E_mat * randn(size(EJ_curr)));
    end
    if apply_noise && isfield(cfg, 'sigma_m_mat') && cfg.sigma_m_mat > 0
        m_curr = m_curr .* (1 + cfg.sigma_m_mat * randn(size(m_curr)));
    end

    %% 2. Applicazione del Danno Strutturale (se presente nello scenario)
    if damage_info.is_damaged && ~isempty(damage_info.element)
        el  = damage_info.element;
        sev = damage_info.severity;
        EJ_curr(el) = EJ_curr(el) .* (1 - sev(:)');
    end

    %% 3. Assemblaggio matrici globali ed estrazione sottomatrici libere (FF)
    [M, K, C] = assem(struct_model.incidenze, struct_model.l, m_curr, ...
                      EA_curr, EJ_curr, struct_model.gamma, struct_model.idb, ...
                      struct_model.conc_data, cfg.rayleigh_ab);
                  
    ndof = struct_model.ndof;
    MFF  = M(1:ndof, 1:ndof);
    KFF  = K(1:ndof, 1:ndof);
    CFF  = C(1:ndof, 1:ndof);

    %% 4. Costruzione posizioni degli assi (Convoglio realistico D1=8m, D2=3m)
    pos_P_in = [];
    x_pos = 0;
    for iv = 1:cfg.n_vagoni
        % Asse anteriore e asse posteriore del vagone corrente
        pos_P_in = [pos_P_in; x_pos; x_pos - cfg.D1];
        % Arretramento al vagone successivo
        x_pos = x_pos - cfg.D1 - cfg.D2;
    end
    
    nP    = length(pos_P_in);
    P_vec = P_val * ones(nP, 1);

    %% 5. Impostazione forzante dinamica e integrazione temporale di Newmark
    if isempty(struct_model.conc_data.deck)
        deck_elements = 1:struct_model.nbeam;
    else
        deck_elements = struct_model.conc_data.deck;
    end

    F_fun = @(tt) Forze_carichi_viaggianti( ...
        tt, ndof, struct_model.xy, struct_model.incidenze, struct_model.l, ...
        struct_model.gamma, struct_model.posiz, deck_elements, pos_P_in, V, P_vec);

    t_sim = (0:cfg.dt_newmark:cfg.t_total)';
    q0      = zeros(ndof, 1);
    q_dot0  = zeros(ndof, 1);
    beta_N  = 1/4; % accelerazione media costante
    gamma_N = 1/2;
    verbose = false; % esecuzione silente in batch

    [~, ~, q_ddot] = Newmark_accelerazioni( ...
        MFF, CFF, KFF, F_fun, t_sim, q0, q_dot0, beta_N, gamma_N, verbose);

    %% 6. Ricampionamento dei 9 sensori alla frequenza target (fs = 200 Hz)
    % Vettore temporale target: esattamente 3000 punti (15 s @ 200 Hz)
    t_target = (0:cfg.pts_total-1)' / cfg.fs_sensor;
    
    raw_sensor_acc = q_ddot(struct_model.sensor_dofs, :); % [9 x nt_sim]
    
    % Interpolazione lineare sul passo target di 200 Hz
    acc_sensors = interp1(t_sim, raw_sensor_acc', t_target, 'linear')'; % [9 x 3000]

    %% 7. Aggiunta rumore di misura del sensore (AWGN)
    if apply_noise && isfield(cfg, 'noise_sensor') && cfg.noise_sensor > 0
        for is = 1:size(acc_sensors, 1)
            sig_std = std(acc_sensors(is, :));
            noise   = cfg.noise_sensor * sig_std * randn(1, cfg.pts_total);
            acc_sensors(is, :) = acc_sensors(is, :) + noise;
        end
    end
end
