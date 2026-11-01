%% RUN_DATA_GENERATION Orchestratore Batch per la Generazione Dati SHM
% Progetto: Master Thesis in Bridge Structural Health Monitoring (BSHM)
% Modello 3: Trave appoggiata da 30m, 9 sensori virtuali, convoglio a vagoni (8m+3m).
%
% Genera due file .mat in output nella cartella data/:
%   1. raw_healthy_runs.mat  -> 154 simulazioni sane (11 carichi x 14 velocità)
%   2. raw_damaged_runs.mat  -> 375 simulazioni con danno (5 posizioni x 5 severità x 3 carichi x 5 velocità)
%
% I dati sono salvati continui a 200 Hz [N_corse x 9 sensori x 3000 punti].
% Il windowing (0.5s) e la normalizzazione verranno effettuati a valle in Python.

clc;
clear;
close all;

%% 1. Inizializzazione percorsi e configurazione
script_dir = fileparts(mfilename('fullpath'));
if ~isempty(script_dir)
    cd(script_dir);
    addpath(script_dir);
end

cfg = config_benchmark();
fprintf('=======================================================\n');
fprintf('  AVVIO GENERAZIONE DATASET SHM - MODELLO 3\n');
fprintf('  Ponte: %s (Luce = %.1f m, 9 Sensori @ %d Hz)\n', cfg.inp_file, cfg.L_bridge, cfg.fs_sensor);
fprintf('=======================================================\n\n');

%% 2. Caricamento della topologia strutturale (UNA SOLA VOLTA)
fprintf('Caricamento struttura da %s...\n', cfg.inp_file);
[file_i, xy, nnod, sizee, idb, ndof, incidenze, l, gamma, ...
 m, EA, EJ, posiz, nbeam, pr, conc_data, rayleigh_ab] = loadstructure(cfg.inp_file);

% Mappatura dei nodi sensore sui corrispondenti GDL verticali (uy)
sensor_dofs = zeros(cfg.n_sensors, 1);
for is = 1:cfg.n_sensors
    node_id = cfg.sensor_nodes(is);
    dof_w = idb(node_id, 2); % 2 = spostamento verticale y
    if dof_w == 0 || dof_w > ndof
        error('Errore: Il nodo sensore %d non ha un GDL verticale valido!', node_id);
    end
    sensor_dofs(is) = dof_w;
end

% Pacchetto dati della struttura
struct_model.incidenze   = incidenze;
struct_model.l           = l;
struct_model.gamma       = gamma;
struct_model.idb         = idb;
struct_model.ndof        = ndof;
struct_model.xy          = xy;
struct_model.posiz       = posiz;
struct_model.nbeam       = nbeam;
struct_model.conc_data   = conc_data;
struct_model.sensor_dofs = sensor_dofs;
struct_model.m_nom       = m;
struct_model.EA_nom      = EA;
struct_model.EJ_nom      = EJ;

fprintf('Topologia verificata con successo: %d nodi, %d travi, %d GDL liberi.\n', nnod, nbeam, ndof);
fprintf('Sensori virtuali (Nodi): %s\n', mat2str(cfg.sensor_nodes));
fprintf('Sensori verticali (GDL): %s\n\n', mat2str(sensor_dofs'));

%% ------------------------------------------------------------------------
%% FASE 1: GENERAZIONE DATASET BASELINE SANO (Training & Validation)
%% ------------------------------------------------------------------------
n_loads  = length(cfg.loads);
n_speeds = length(cfg.speeds_ms);
n_healthy_runs = n_loads * n_speeds;

fprintf('--- 1. Generazione Baseline Sana: %d corse totali ---\n', n_healthy_runs);
fprintf('    (%d carichi [%.0f - %.0f kN] x %d velocità [%.0f - %.0f km/h])\n', ...
    n_loads, min(cfg.loads)/1e3, max(cfg.loads)/1e3, ...
    n_speeds, min(cfg.speeds_kmh), max(cfg.speeds_kmh));

acc_healthy  = zeros(n_healthy_runs, cfg.n_sensors, cfg.pts_total, 'single');
meta_healthy = zeros(n_healthy_runs, 7); 
% Colonne meta: [RunID, Load_N, Speed_kmh, Speed_ms, IsDamaged (0), DamageElement (0), DamageSeverity (0)]

t_start_healthy = tic;
run_idx = 1;

for iL = 1:n_loads
    P_val = cfg.loads(iL);
    for iV = 1:n_speeds
        V_ms  = cfg.speeds_ms(iV);
        V_kmh = cfg.speeds_kmh(iV);
        
        dmg_info.is_damaged = false;
        dmg_info.element    = [];
        dmg_info.severity   = 0;
        
        % Simulazione singolo transito con rumore sensori (5%) e materiale (2%)
        [acc_run, t_target] = simulate_transit(cfg, struct_model, V_ms, P_val, dmg_info, true);
        
        acc_healthy(run_idx, :, :) = single(acc_run);
        meta_healthy(run_idx, :)   = [run_idx, P_val, V_kmh, V_ms, 0, 0, 0];
        
        if mod(run_idx, 20) == 0 || run_idx == n_healthy_runs
            fprintf('  [Sano] Progresso: %d/%d corse completate (%.1f%%)...\n', ...
                run_idx, n_healthy_runs, 100 * run_idx / n_healthy_runs);
        end
        run_idx = run_idx + 1;
    end
end
time_healthy = toc(t_start_healthy);
fprintf('✓ Generazione Baseline Sana completata in %.2f s (media: %.3f s/corsa).\n', ...
    time_healthy, time_healthy / n_healthy_runs);

% Salvataggio Baseline Sana su file .mat
file_out_healthy = fullfile(cfg.output_dir, 'raw_healthy_runs.mat');
save(file_out_healthy, 'acc_healthy', 'meta_healthy', 't_target', 'cfg', '-v7.3');
fprintf('  File salvato: %s\n\n', file_out_healthy);

%% ------------------------------------------------------------------------
%% FASE 2: GENERAZIONE DATASET DANNEGGIATO (Test Set per Anomaly Detection)
%% ------------------------------------------------------------------------
n_dmg_elements   = length(cfg.damage_elements);
n_dmg_severities = length(cfg.damage_severities);
n_dmg_loads      = length(cfg.damage_loads);
n_dmg_speeds     = length(cfg.damage_speeds_ms);

n_damaged_runs = n_dmg_elements * n_dmg_severities * n_dmg_loads * n_dmg_speeds;

fprintf('--- 2. Generazione Scenari di Danno: %d corse totali ---\n', n_damaged_runs);
fprintf('    (%d elementi x %d severità x %d carichi x %d velocità)\n', ...
    n_dmg_elements, n_dmg_severities, n_dmg_loads, n_dmg_speeds);

acc_damaged  = zeros(n_damaged_runs, cfg.n_sensors, cfg.pts_total, 'single');
meta_damaged = zeros(n_damaged_runs, 7);

t_start_damaged = tic;
run_idx = 1;

for i_el = 1:n_dmg_elements
    el_id = cfg.damage_elements(i_el);
    for i_sev = 1:n_dmg_severities
        sev = cfg.damage_severities(i_sev);
        
        dmg_info.is_damaged = true;
        dmg_info.element    = el_id;
        dmg_info.severity   = sev;
        
        for iL = 1:n_dmg_loads
            P_val = cfg.damage_loads(iL);
            for iV = 1:n_dmg_speeds
                V_ms  = cfg.damage_speeds_ms(iV);
                V_kmh = cfg.damage_speeds_kmh(iV);
                
                [acc_run, ~] = simulate_transit(cfg, struct_model, V_ms, P_val, dmg_info, true);
                
                acc_damaged(run_idx, :, :) = single(acc_run);
                meta_damaged(run_idx, :)   = [run_idx, P_val, V_kmh, V_ms, 1, el_id, sev];
                
                if mod(run_idx, 50) == 0 || run_idx == n_damaged_runs
                    fprintf('  [Danno] Progresso: %d/%d corse completate (%.1f%%)...\n', ...
                        run_idx, n_damaged_runs, 100 * run_idx / n_damaged_runs);
                end
                run_idx = run_idx + 1;
            end
        end
    end
end
time_damaged = toc(t_start_damaged);
fprintf('✓ Generazione Scenari di Danno completata in %.2f s (media: %.3f s/corsa).\n', ...
    time_damaged, time_damaged / n_damaged_runs);

% Salvataggio Dataset Danneggiato su file .mat
file_out_damaged = fullfile(cfg.output_dir, 'raw_damaged_runs.mat');
save(file_out_damaged, 'acc_damaged', 'meta_damaged', 't_target', 'cfg', '-v7.3');
fprintf('  File salvato: %s\n\n', file_out_damaged);

fprintf('=======================================================\n');
fprintf('  TUTTE LE SIMULAZIONI COMPLETATE CON SUCCESSO!\n');
fprintf('  1. Dati Sani:       %s [%d x 9 x 3000]\n', file_out_healthy, n_healthy_runs);
fprintf('  2. Dati Danneggiati: %s [%d x 9 x 3000]\n', file_out_damaged, n_damaged_runs);
fprintf('=======================================================\n');
