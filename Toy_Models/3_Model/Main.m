%% Unified 2D Beam/Truss Railway Bridge Dynamic Simulation - Model 3
% Progetto: Master Thesis in Bridge Structural Health Monitoring (BSHM)
% Modello unificato 2D: elementi trave Euler-Bernoulli + assiale, elementi concentrati,
% carichi da convoglio ferroviario ad assi multipli ed integrazione temporale Newmark.

clc;
close all;

%% 1. Caricamento della struttura dal file .inp
if ~exist('nome_file_inp', 'var') || isempty(nome_file_inp)
    nome_file_inp = 'bridge_beam.inp'; % File di input di default
end

fprintf('=== Caricamento modello da file: %s ===\n', nome_file_inp);
[file_i, xy, nnod, sizee, idb, ndof, incidenze, l, gamma, m, EA, EJ, posiz, nbeam, pr, conc_data] = ...
    loadstructure(nome_file_inp);

fprintf('Nodi totali: %d | Travi: %d | GDL liberi: %d\n', nnod, nbeam, ndof);

%% 2. Visualizzazione della struttura indeformata e condizioni di vincolo
dis_stru(posiz, l, gamma, xy, pr, idb, ndof);

%% 3. Assemblaggio matrici globali ed estrazione sottomatrici libere (FF)
% Parametri smorzamento di Rayleigh: C = alpha * M + beta * K + C_conc
rayleigh_ab = [0.1, 2e-4];

[M, K, C] = assem(incidenze, l, m, EA, EJ, gamma, idb, conc_data, rayleigh_ab);

% Partizionamento sui soli gradi di libertà liberi
MFF = M(1:ndof, 1:ndof);
KFF = K(1:ndof, 1:ndof);
CFF = C(1:ndof, 1:ndof);

%% 4. Analisi modale
n_modi_interesse = min(7, ndof);

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
n_plot_modi = min(3, n_modi_interesse);
modes_plot_conc_element(phi_w, xy, freq, n_plot_modi, [(1:nnod-1)', (2:nnod)']);

% Plot della configurazione deformata 2D del primo modo
figure('Name', 'Forma Modale 1 (Configurazione 2D)', 'NumberTitle', 'off');
diseg2(Phi(:, 1), 1.5, incidenze, l, gamma, posiz, idb, xy);
title(sprintf('Modo 1 - Frequenza: %.3f Hz', freq(1)));

%% 6. Definizione del convoglio ferroviario (Moving Train)
% Identificazione elementi appartenenti al piano del ferro (Deck)
if isempty(conc_data.deck)
    deck_elements = 1:nbeam; % se non specificato, si assume tutta la travata
else
    deck_elements = conc_data.deck;
end

V       = 40.0;                       % Velocità convoglio [m/s] (~144 km/h)
n_P     = 10;                         % Numero di assi del convoglio
d_P     = 5.0;                        % Interasse tra assi consecutivi [m]
P_val   = 1.2e5;                      % Carico per asse [N] (~12 tonnellate)

P_vec    = P_val * ones(n_P, 1);
pos_P_in = -(0:n_P-1)' * d_P;         % Posizioni iniziali a monte del ponte [m]

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

dt = 1e-4;                                 % Passo temporale di integrazione [s]
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
    title(sprintf('Accelerazione Verticale Nodo %d (Sensore SHM)', node_mid));
end

fprintf('\n=== Simulazione completata con successo! ===\n');
