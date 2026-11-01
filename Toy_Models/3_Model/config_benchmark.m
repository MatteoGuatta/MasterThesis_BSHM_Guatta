function cfg = config_benchmark()
% CONFIG_BENCHMARK Parametri e configurazione del benchmark Beam-Bridge
% per l'addestramento e il test dell'Anomaly Detection (Master Thesis SHM)

    %% 1. File del modello FEM
    cfg.inp_file = 'bridge_beam.inp';
    
    %% 2. Sensori Virtuali
    % 9 sensori equispaziati a 3m: da x = 3m (nodo 4) a x = 27m (nodo 28)
    cfg.sensor_nodes = [4, 7, 10, 13, 16, 19, 22, 25, 28];
    cfg.n_sensors    = length(cfg.sensor_nodes);
    
    %% 3. Geometria del Convoglio Ferroviario (Vagoni 8m + 3m)
    cfg.n_vagoni = 8;     % Numero di vagoni
    cfg.D1       = 8.0;   % Interasse interno tra i 2 assi del medesimo vagone [m]
    cfg.D2       = 3.0;   % Distanza tra asse posteriore e asse anteriore tra vagoni limitrofi [m]
    
    %% 4. Griglie EOV (Variabilità Operativa)
    % 11 classi di carico per asse: da 80 kN a 130 kN (step 5 kN)
    cfg.loads = (80:5:130) * 1e3; % [N]
    
    % 14 classi di velocità: da 38 km/h a 77 km/h (step 3 km/h)
    cfg.speeds_kmh = (38:3:77);          % [km/h]
    cfg.speeds_ms  = cfg.speeds_kmh / 3.6; % [m/s]
    
    %% 5. Parametri Temporali e Campionamento
    cfg.fs_sensor  = 200;                                % Frequenza target sensori [Hz]
    cfg.t_total    = 15.0;                               % Durata totale del transito [s]
    cfg.pts_total  = round(cfg.t_total * cfg.fs_sensor); % 3000 punti per sensore
    cfg.dt_newmark = 5e-4;                               % Passo di integrazione Newmark [s] (0.5 ms: 4x più veloce, errore < 0.25%)
    
    %% 6. Parametri di Smorzamento di Rayleigh [alpha, beta]
    cfg.rayleigh_ab = [0.1, 2e-4];
    
    %% 7. Incertezze e Rumore
    cfg.noise_sensor = 0.05; % 5% rumore bianco gaussiano additivo (AWGN)
    cfg.sigma_E_mat  = 0.02; % 2% variabilità naturale del modulo elastico tra elementi
    cfg.sigma_m_mat  = 0.02; % 2% variabilità naturale della massa distribuita
    
    %% 8. Scenari di Danno per il Test Set
    % Posizioni spaziali del danno (elementi trave da 1m lungo i 30m di luce)
    cfg.damage_elements   = [5, 10, 15, 20, 25];         % x = 5m, 10m, 15m (mezzeria), 20m, 25m
    cfg.damage_severities = [0.05, 0.10, 0.25, 0.50, 0.75]; % Riduzione di EJ (-5%, -10%, -25%, -50%, -75%)
    
    % Sottoinsieme di condizioni operative per i test di danno
    cfg.damage_speeds_kmh = [38, 44, 51, 62, 77];        % velocità rappresentative (risonanze e valli)
    cfg.damage_speeds_ms  = cfg.damage_speeds_kmh / 3.6;
    cfg.damage_loads      = [80, 105, 130] * 1e3;        % carichi: min, medio, max [N]
    
    %% 9. Cartella di output dei dati generati
    script_dir = fileparts(mfilename('fullpath'));
    cfg.output_dir = fullfile(script_dir, 'data');
    if ~exist(cfg.output_dir, 'dir')
        mkdir(cfg.output_dir);
    end
end
