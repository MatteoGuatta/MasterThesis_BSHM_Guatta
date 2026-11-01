function cfg = config_benchmark()
% CONFIG_BENCHMARK Raccoglie tutti i parametri operativi, geometrici 
% e di simulazione per il benchmark del ponte ferroviario.

    %% 1. File di Input Radice
    % È il punto di partenza: tutto il modello geometrico e le proprietà 
    % della trave verranno letti da questo file tramite loadstructure()
    cfg.inp_file = 'bridge_beam.inp';

    [file_i, xy, nnod, sizee, idb, ndof, incidenze, l, gamma, m, EA, EJ, posiz, nbeam, pr, conc_data] = loadstructure(cfg.inp_file);

    cfg.xy = xy;
    cfg.nnod = nnod;
    cfg.sizee = sizee;
    cfg.idb = idb;
    cfg.ndof = ndof;
    cfg.incidenze = incidenze;
    cfg.l = l;
    cfg.gamma = gamma;
    cfg.m = m;
    cfg.EA = EA;
    cfg.EJ = EJ;
    cfg.posiz = posiz;
    cfg.nbeam = nbeam;
    cfg.pr = pr;
    cfg.conc_data = conc_data;

    %% 2. Sensori Virtuali & Deck elements
    % Definiamo i 9 nodi corrispondenti alle posizioni x = [3, 6, 9, 12, 15, 18, 21, 24, 27] m
    cfg.sensor_nodes = [4, 7, 10, 13, 16, 19, 22, 25, 28];
    cfg.n_sensors    = length(cfg.sensor_nodes);

    if isempty(conc_data.deck)
        cfg.deck_elements = 1:cfg.nbeam; % se non specificato, si assume tutta la travata
    else
        cfg.deck_elements = conc_data.deck;
    end

    %% 3. Geometria del Convoglio (Discretizzazione a vagoni)
    cfg.n_vagoni = 8;     % Numero di carrozze del treno
    cfg.D1       = 8.0;   % Distanza tra i due assi dello stesso vagone [m]
    cfg.D2       = 3.0;   % Distanza tra l'asse posteriore e anteriore tra vagoni consecutivi [m]

    %% 4. Variabilità Operativa (EOV)
    % 11 classi di carico: da 80 kN a 130 kN a passi di 5 kN
    cfg.loads = (80:5:130) * 1e3; % [N]
    
    % 14 classi di velocità: da 38 km/h a 77 km/h a passi di 3 km/h
    cfg.speeds_kmh = (38:3:77);          % [km/h]
    cfg.speeds_ms  = cfg.speeds_kmh / 3.6; % [m/s]

    %% 5. Parametri Temporali e Campionamento
    cfg.fs_sensor  = 200;                                % Frequenza campionamento target dei sensori [Hz]
    cfg.t_total    = 15.0;                               % Durata della registrazione di ogni corsa [s]
    cfg.pts_total  = round(cfg.t_total * cfg.fs_sensor); % Esattamente 3000 punti per sensore
    cfg.dt_newmark = 5e-3;                               % Passo temporale Newmark (0.5 ms per accuratezza e velocità)
    
    %% 6. Smorzamento Strutturale (Rayleigh)
    % C = alpha * M + beta * K
    cfg.rayleigh_ab = [0.1, 2e-4];

    %% 7. Incertezze e Rumore
    cfg.noise_sensor = 0.00; % 5% rumore bianco di misura sui sensori (AWGN)
    cfg.sigma_E_mat  = 0.02; % 2% variabilità naturale del modulo elastico del materiale
    cfg.sigma_m_mat  = 0.02; % 2% variabilità naturale della massa distribuita

    %% 8. Scenari di Danno per il Test Set
    
    cfg.scenari_danno = {
    % --- TEST 1: Danno Locale Incrementale ---
    'Locale_Quarto_L4',   'locale',    8,         [0.10, 0.25, 0.50, 0.75];
    'Locale_Mezzeria',    'locale',    15,        [0.10, 0.25, 0.50, 0.75];
    'Locale_TreQuarti',   'locale',    23,        [0.10, 0.25, 0.50, 0.75];
    
    % --- TEST 2: Danno Globale ---
    'Degrado_Globale',    'globale',   1:30,      0.30; % -20% su tutti gli elementi
    
    % --- TEST 3: Doppio Danno Simultaneo ---
    'Doppio_Difetto',     'multiplo',  [8, 22],   [0.25, 0.25]; % -25% su entrambi
    };
    
    % Condizioni operative ridotte per testare il danno (senza esplodere con troppe corse)
    cfg.damage_speeds_kmh = [38, 56,77];            % Risonanze e valli
    cfg.damage_speeds_ms  = cfg.damage_speeds_kmh / 3.6;
    cfg.damage_loads      = [80, 105, 130] * 1e3;            % Minimo, medio, massimo carico [N]

    %% 9. Cartella di Salvataggio Dati
    script_dir = fileparts(mfilename('fullpath'));
    cfg.output_dir = fullfile(script_dir, 'data');
    if ~exist(cfg.output_dir, 'dir')
        mkdir(cfg.output_dir);
    end
end