clc;
clear;
close all;


%% Lettura del file di configurazione

cfg = config_benchmark();

script_dir = fileparts(mfilename('fullpath'));
if ~isempty(script_dir)
    cd(script_dir);
    addpath(script_dir);
end


% Stampa dati relativi a struttura e sensoristica

fprintf('Topologia verificata con successo: %d nodi, %d travi, %d GDL liberi.\n', cfg.nnod, cfg.nbeam, cfg.ndof);
fprintf('Sensori virtuali (Nodi): %s\n', mat2str(cfg.sensor_nodes));
%fprintf('Sensori verticali (GDL): %s\n\n', mat2str(cfg.sensor_dofs'));


%% 1. Generazione Dataset Baseline Sano

n_loads = length(cfg.loads);
n_speeds = length(cfg.speeds_kmh);
n_healthy_runs = n_speeds * n_loads;

fprintf('--- 1. Generazione Baseline Sana: %d corse totali ---\n', n_healthy_runs);
fprintf('    (%d carichi [%.0f - %.0f kN] x %d velocità [%.0f - %.0f km/h])\n', ...
    n_loads, min(cfg.loads)/1e3, max(cfg.loads)/1e3, ...
    n_speeds, min(cfg.speeds_kmh), max(cfg.speeds_kmh));


acc_healthy = zeros(n_healthy_runs, cfg.n_sensors, cfg.pts_total);
% colonne meta_healthy = [RunID, Load_N, Speed_kmh, Speed_ms, IsDamaged (0), DamageElement (0), DamageSeverity (0)]
meta_healthy = zeros(n_healthy_runs, 7);

run_idx = 1;
t_start_healthy = tic;

for i = 1:n_loads
    P_val = cfg.loads(i);
    for k = 1:n_speeds
        V = cfg.speeds_ms(k); % velocity expressed in [m/s]

        dmg_info.is_damaged = false;
        dmg_info.element = [];
        dmg_info.severity = 0;

        % Simulazione Singolo transito 

        [acc_run, t_target] = simulate_transit(cfg, V, P_val, dmg_info);

        acc_healthy(run_idx,:,:) = acc_run;
        meta_healthy(run_idx,:) = [run_idx, P_val, V*3.6, V, 0, 0, 0];

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


%% 2. Generazione Dataset Danno

% Selezione manuale della casistica di danno da simulare [1-Locale_Quarto_L4 , 2 Locale_Mezzeria, 3 Locale_TreQuarti, 4 Degrado_Globale, 5 Doppio_Difetto]
idx_scenario = 1;

nome_scenario = cfg.scenari_danno{idx_scenario, 1}; % 
tipo_scenario = cfg.scenari_danno{idx_scenario, 2}; % 
el_target     = cfg.scenari_danno{idx_scenario, 3}; % 
severita_vec  = cfg.scenari_danno{idx_scenario, 4}; % 

% Calcolo del numero di simulazioni
n_speeds_damage = length(cfg.damage_speeds_kmh);
n_loads_damage = length(cfg.damage_loads);
n_severities = length(severita_vec);

n_damaged_runs = n_speeds_damage * n_loads_damage * n_severities;


fprintf('--- 1. Generazione Dataset Danno: %d corse totali ---\n', n_damaged_runs);
fprintf('    (%d carichi [%.0f - %.0f kN] x %d velocità [%.0f - %.0f km/h])\n', ...
    n_loads_damage, min(cfg.damage_loads)/1e3, max(cfg.damage_loads)/1e3, ...
    n_speeds_damage, min(cfg.damage_speeds_kmh), max(cfg.damage_speeds_kmh));

acc_damaged = zeros(n_damaged_runs, cfg.n_sensors, cfg.pts_total);
meta_damaged = zeros(n_damaged_runs, 7);

% setup simulazione
run_idx = 1;
t_start_damaged = tic;

for ii = 1:n_severities
    for i = 1:n_loads_damage
        P_val = cfg.damage_loads (i);
        for k = 1:n_speeds_damage
            V = cfg.damage_speeds_ms(k); % velocity expressed in [m/s]
    
            dmg_info.is_damaged = true;
            dmg_info.element = el_target;
            dmg_info.severity = severita_vec(ii);
    
            % Simulazione Singolo transito 
    
            [acc_run, t_target] = simulate_transit(cfg, V, P_val, dmg_info);
    
            acc_damaged(run_idx,:,:) = acc_run;
            meta_damaged(run_idx, :) = [run_idx, P_val, V*3.6, V, 1, idx_scenario, severita_vec(ii)];    
            
            if mod(run_idx, 20) == 0 || run_idx == n_damaged_runs
                fprintf('  [Damaged] Progresso: %d/%d corse completate (%.1f%%)...\n', ...
                    run_idx, n_damaged_runs, 100 * run_idx / n_damaged_runs);
            end
            run_idx = run_idx + 1;
    
        end
    end
end

time_damaged = toc(t_start_damaged);
fprintf('✓ Generazione dataset danno in %.2f s (media: %.3f s/corsa).\n', ...
    time_damaged, time_damaged / n_damaged_runs);

% Salvataggio Baseline Sana su file .mat
file_out_damaged = fullfile(cfg.output_dir, 'raw_damaged_runs.mat');
save(file_out_damaged, 'acc_damaged', 'meta_damaged', 't_target', 'cfg', '-v7.3');
fprintf('  File salvato: %s\n\n', file_out_damaged);

