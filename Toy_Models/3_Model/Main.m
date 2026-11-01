%% Unified 2D Beam/Truss Railway Bridge Dynamic Simulation - Model 3
% Progetto: Master Thesis in Bridge Structural Health Monitoring (BSHM)
% Modello unificato 2D: elementi trave Euler-Bernoulli + assiale, elementi concentrati,
% carichi da convoglio ferroviario ad assi multipli ed integrazione temporale Newmark.

clc;
close all;

% Aggiunge la cartella dello script al path di MATLAB
cartella_script = fileparts(mfilename('fullpath'));
if ~isempty(cartella_script)
    addpath(cartella_script);
end

%% 1. Caricamento della struttura baseline dal file .inp
fprintf('=== Caricamento modello da file: %s ===\n', 'bridge_beam.inp');
[file_i, xy, nnod, sizee, idb, ndof, incidenze, l, gamma, m, EA, EJ, posiz, nbeam, pr, conc_data] = ...
    loadstructure('bridge_beam.inp');

fprintf('Nodi totali: %d | Travi: %d | GDL liberi: %d\n', nnod, nbeam, ndof);

%% 2. Visualizzazione della struttura indeformata e condizioni di vincolo
dis_stru(posiz, l, gamma, xy, pr, idb, ndof);

%% Definizione Scenari di Danno
scenarios = {
    % Nome scenario,        Elementi target,      Fattori di riduzione
    'Baseline_Sano',         [],                  [];
    %'Danno_Locale_Singolo',  13,                  0.50;
    %'Danno_Multiplo',        [13, 14],            [0.50, 0.50];
    %'Danno_Diffuso_Campata', 15:20,               repmat(0.30, 1, 6);
    %'Corrosione_Globale',    1:length(l),         repmat(0.25, 1, length(l));
};

%% Loop per generazione e simulazione batch
for s = 1:size(scenarios, 1)
    scen_name = scenarios{s, 1};
    el_id     = scenarios{s, 2}; % Vettore indici elementi
    delta     = scenarios{s, 3}; % Vettore percentuali riduzione (es. 0.20)
    
    % Reset alle proprietà nominali intatte
    EJ_curr = EJ;
    EA_curr = EA;
    
    % Applicazione del danno (se presente nello scenario)
    if ~isempty(el_id)
        EJ_curr(el_id) = EJ_curr(el_id) .* (1 - delta(:)');
        EA_curr(el_id) = EA_curr(el_id) .* (1 - delta(:)');
    end
   
    %% 3. Assemblaggio matrici globali ed estrazione sottomatrici libere (FF)
    % Parametri smorzamento di Rayleigh: C = alpha * M + beta * K + C_conc
    rayleigh_ab = [0.1, 2e-4];
    
    % B = [0.01 0.015 0.018]';
    % A = zeros(3,2);
    % for ii=1:3
    % A(ii,:) = [1/(2*omega(ii)) omega(ii)/2];
    % end
    % ab = A\B;
    
    [M, K, C] = assem(incidenze, l, m, EA_curr, EJ_curr, gamma, idb, conc_data, rayleigh_ab);
    
    % Partizionamento sui soli gradi di libertà liberi
    MFF = M(1:ndof, 1:ndof);
    KFF = K(1:ndof, 1:ndof);
    CFF = C(1:ndof, 1:ndof);
    
    %% 4. Analisi modale
    n_modi_interesse = min(8, ndof);
    
    [freq, Phi, zeta, tab_modi] = modal_analysis(MFF, KFF, CFF, n_modi_interesse);
    
    fprintf('\n--- Risultati Analisi Modale ---\n');
    disp(tab_modi);
    
    %% 5. Plot delle forme modali
    % Estrazione componenti verticali per il plot 1D
    phi_w = zeros(nnod, n_modi_interesse);
    for in = 1:nnod
        dof_w = idb(in, 2); % Grado di libertà verticale (y)
        if dof_w > 0 && dof_w <= ndof
            phi_w(in, :) = Phi(dof_w, 1:n_modi_interesse);
        end
    end
    
    % Plot profilo verticale dei modi
    n_plot_modi = min(7, n_modi_interesse);
    modes_plot_conc_element(phi_w, xy, freq, n_plot_modi, [(1:nnod-1)', (2:nnod)']);
    
    % Plot della configurazione deformata 2D del primo modo
    figure('Name', 'Forma Modale 1 (Configurazione 2D)', 'NumberTitle', 'off');
    diseg2(Phi(:, 1), 1.5, incidenze, l, gamma, posiz, idb, xy);
    title(sprintf('Modo 1 - Frequenza: %.3f Hz', freq(1)));
    
    % plot per evidenziare la prima frequenza assiale 
    figure('Name', 'Forma Modale 1 (Configurazione 2D)', 'NumberTitle', 'off');
    diseg2(Phi(:, 2), 1.5, incidenze, l, gamma, posiz, idb, xy);
    title(sprintf('Modo 1 - Frequenza: %.3f Hz', freq(2)));

    
    %% 6. Definizione del convoglio ferroviario (Moving Train)
    % Identificazione elementi appartenenti al piano del ferro (Deck)
    if isempty(conc_data.deck)
        deck_elements = 1:nbeam; % se non specificato, si assume tutta la travata
    else
        deck_elements = conc_data.deck;
    end
    
    V = 10.55;       % velocita' [m/s]
    n_vagoni = 15;
    n_P = 30;                % numero di carichi
    d_P = 8;      % Interasse tra ruote del medesimo vagone [m]
    d_ext =3;     % Interasse tra asse anteriore e posteriore di vagoni limitrofi
    P_val = 1e5;  % intensita' di ogni carico [N]
    
    P_vec    = P_val * ones(n_P, 1);
    pos_P_in = [];
    x_curr = 0;

    for iv = 1:n_vagoni
        % Asse 1 e Asse 2 del vagone iv
        pos_P_in = [pos_P_in; -x_curr; -(x_curr + d_P)];
        
        % Avanzamento al vagone successivo
        x_curr = x_curr + d_ext + d_P;
    end
    
    %% 7. Intervallo temporale di integrazione
    deck_el      = deck_elements(:);
    x_start_deck = posiz(deck_el, 1);
    l_deck       = l(deck_el); l_deck = l_deck(:);
    gamma_deck   = gamma(deck_el); gamma_deck = gamma_deck(:);
    x_end_deck   = x_start_deck + l_deck .* cos(gamma_deck);
    
    x_deck_start = min([x_start_deck; x_end_deck]);
    x_deck_end   = max([x_start_deck; x_end_deck]);
    
    t_exit = (x_deck_end - min(pos_P_in)) / V; % Istante di uscita dell'ultimo asse
    t_free = 1.5;                              % Secondi di vibrazione libera successivi
    t_end  = t_exit + t_free;
    
    dt = 5e-3;                                 % Passo temporale di integrazione [s] -- (500 Hz) / inverso della frequenza di campionamento dei sensori
    n_steps = ceil(t_end / dt);
    t = linspace(0, t_end, n_steps + 1)';
    
    fprintf('\n=== Avvio Simulazione Dinamica Passaggio Treno ===\n');
    fprintf('Durata transito: %.2f s | Tempo totale: %.2f s | Timestep dt: %.1e s\n', t_exit, t_end, dt);
    
    % Funzione forzante dinamica
    F_fun = @(tt) Forze_carichi_viaggianti( ...
        tt, ndof, xy, incidenze, l, gamma, posiz, deck_elements, pos_P_in, V, P_vec);
    
    %% 8. Integrazione diretta con algoritmo di Newmark
    beta_N  = 1/4; % Variante accelerazione media costante (incondizionatamente stabile)
    gamma_N = 1/2;
    
    q0     = zeros(ndof, 1);
    q_dot0 = zeros(ndof, 1);
    
    [q, q_dot, q_ddot] = Newmark_accelerazioni( ...
        MFF, CFF, KFF, F_fun, t, q0, q_dot0, beta_N, gamma_N);
    
    %% 9. Estrazione e Plot Risposta Dinamica (Sensore a mezzeria)
    % Individuazione del nodo più vicino a mezzeria dell'impalcato
    x_mid = 0.5 * (x_deck_start + x_deck_end);
    [~, node_mid] = min(abs(xy(:, 1) - x_mid));
    dof_mid_y = idb(node_mid, 2); % GDL verticale (uy)
   

    % Nodi fisici su cui sono posizionati i sensori verticali
    
    sensors_node = 4:3:28;
    n_sensors = length(sensors_node);

    % Estrazione dei gradi di libertà verticali (colonna 2 = traslazione y)
    sensors_dof = zeros(1,n_sensors);
    sensors_acc = zeros(n_sensors, n_steps+1);

    for i = 1:length(sensors_node)
        sensors_dof(i) = idb(sensors_node(i), 2);
        sensors_acc(i,:) = q_ddot(sensors_dof(i), :);
    end

    if dof_mid_y > 0 && dof_mid_y <= ndof
        disp_mid = q(dof_mid_y, :)';       % [m]
        acc_mid  = q_ddot(dof_mid_y, :)';  % [m/s^2]
    
        figure('Name', 'Risposta Dinamica a Mezzeria', 'NumberTitle', 'off');
    
        % Spostamento
        subplot(2, 1, 1);
        plot(t, disp_mid * 1000, 'b', 'LineWidth', 1.5);
        grid on;
        xlabel('Tempo [s]');
        ylabel('Freccia u_y [mm]');
        title(sprintf('Spostamento Verticale Nodo %d (Mezzeria X = %.2f m)', node_mid, xy(node_mid, 1)));
    
        % Accelerazione (segnale sensore SHM)
        subplot(2, 1, 2);
        plot(t, acc_mid, 'r', 'LineWidth', 1.2);
        grid on;
        xlabel('Tempo [s]');
        ylabel('$\ddot{u}_y \ [\mathrm{m/s^2}]$', 'Interpreter', 'latex');
        title(sprintf('Accelerazione Verticale Nodo %d', node_mid));
    end
    
    for i = 1:length(sensors_node)
        figure; hold on
        plot(t, sensors_acc(i,:), 'b', 'LineWidth', 1.5);
        grid on;
        xlabel('Tempo [s]');
        ylabel('$\ddot{u}_y \ [\mathrm{m/s^2}]$', 'Interpreter', 'latex');
        title(sprintf('Accelerazione Sensore %d (X = %.2f m)', sensors_node(i), xy(sensors_node(i), 1)))
    end
    fprintf('\n=== Simulazione completata con successo! ===\n');

end