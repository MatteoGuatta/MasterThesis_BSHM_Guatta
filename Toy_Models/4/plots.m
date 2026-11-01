clc;
clear;
close all;

%% 1. Caricamento dati
scriptPath = fileparts(mfilename('fullpath'));

fileHealthy = fullfile(scriptPath, 'data', 'raw_healthy_runs.mat');
fileDamaged = fullfile(scriptPath, 'data', 'raw_damaged_runs_global.mat');

healthy = load(fileHealthy);
damaged = load(fileDamaged);

acc_healthy  = healthy.acc_healthy;
meta_healthy = healthy.meta_healthy;

acc_dmg  = damaged.acc_damaged;
meta_dmg = damaged.meta_damaged;

t_healthy = healthy.t_target(:);
t_dmg     = damaged.t_target(:);

%% 2. Parametri dei confronti
sensor = 3; % Nodo 10, x = 9 m

v = [38, 56, 77];             % [km/h]
P = [80, 105, 130] * 1e3;     % Carico per asse [N]
severity = 0.30;

V_fixed = 56;                % [km/h]
P_fixed = 105e3;             % [N]
scenario = 4;               % Danno sull'elemento 8

% Scegli: "first", "velocità", "carico", "danno"
index_plot = "danno";

%% 3. Creazione figura
switch index_plot

    case "first"
        % Transito sano di riferimento
        idx = find_run(meta_healthy, P_fixed, V_fixed);

        signal = squeeze(acc_healthy(idx, sensor, :));

        figure('Color', 'w');
        plot(t_healthy, signal, 'LineWidth', 1);

        xlabel('Tempo [s]');
        ylabel('Accelerazione [m/s^2]');
        title(sprintf( ...
            'Sano — nodo 10 — %.0f kN/asse — %.0f km/h', ...
            P_fixed / 1e3, V_fixed));
        grid on;
        xlim([0, healthy.cfg.t_total]);

        % Uscita dell'ultimo asse dal ponte
        trainLength = (healthy.cfg.n_vagoni - 1) * ...
            (healthy.cfg.D1 + healthy.cfg.D2) + healthy.cfg.D1;

        deck = healthy.cfg.deck_elements;
        deck = deck(:);

        xStart = healthy.cfg.posiz(deck, 1);
        lengths = healthy.cfg.l(deck);
        angles = healthy.cfg.gamma(deck);

        xEnd = xStart + lengths(:) .* cos(angles(:));
        xDeckEnd = max([xStart; xEnd]);

        tExit = (xDeckEnd + trainLength) / (V_fixed / 3.6);

        xline(tExit, '--k', 'Uscita ultimo asse', ...
            'LabelVerticalAlignment', 'bottom');

    case "velocità"
        % Struttura sana, carico fisso, velocità variabile
        figure('Color', 'w');
        layout = tiledlayout(numel(v), 1, ...
            'TileSpacing', 'compact', 'Padding', 'compact');

        ax = gobjects(numel(v), 1);
        maxAmplitude = 0;

        for k = 1:numel(v)
            idx = find_run(meta_healthy, P_fixed, v(k));
            signal = squeeze(acc_healthy(idx, sensor, :));

            ax(k) = nexttile(layout);
            plot(t_healthy, signal, 'LineWidth', 1);

            title(sprintf('Velocità: %.0f km/h', v(k)));
            ylabel('a [m/s^2]');
            grid on;

            maxAmplitude = max(maxAmplitude, max(abs(signal)));
        end

        linkaxes(ax, 'xy');
        xlim(ax(1), [0, healthy.cfg.t_total]);

        if maxAmplitude > 0
            ylim(ax(1), 1.05 * [-maxAmplitude, maxAmplitude]);
        end

        xlabel(ax(end), 'Tempo [s]');
        title(layout, sprintf( ...
            'Effetto della velocità — sano — nodo 10 — %.0f kN/asse', ...
            P_fixed / 1e3));

       case "carico"
        % Struttura sana, velocità fissa, carico variabile
        figure('Color', 'w');
        layout = tiledlayout(numel(P), 1, ...
            'TileSpacing', 'compact', 'Padding', 'compact');

        ax = gobjects(numel(P), 1);
        colors = lines(numel(P));
        maxAmplitude = 0;

        for k = 1:numel(P)
            idx = find_run(meta_healthy, P(k), V_fixed);
            signal = squeeze(acc_healthy(idx, sensor, :));

            ax(k) = nexttile(layout);
            plot(t_healthy, signal, ...
                'Color', colors(k, :), ...
                'LineWidth', 1);

            title(sprintf('Carico: %.0f kN/asse', P(k) / 1e3));
            ylabel('a [m/s^2]');
            grid on;

            maxAmplitude = max(maxAmplitude, max(abs(signal)));
        end

        linkaxes(ax, 'xy');
        xlim(ax(1), [0, healthy.cfg.t_total]);

        if maxAmplitude > 0
            ylim(ax(1), 1.05 * [-maxAmplitude, maxAmplitude]);
        end

        xlabel(ax(end), 'Tempo [s]');
        title(layout, sprintf( ...
            'Effetto del carico — sano — nodo 10 — %.0f km/h', ...
            V_fixed));

    case "danno"
        % Carico e velocità fissi, severità variabile
        % Primo pannello: sano. Successivi: danneggiati.
        figure('Color', 'w');
        layout = tiledlayout(numel(severity) + 1, 1, ...
            'TileSpacing', 'compact', 'Padding', 'compact');

        ax = gobjects(numel(severity) + 1, 1);
        colors = lines(numel(severity));
        maxAmplitude = 0;

        idx = find_run(meta_healthy, P_fixed, V_fixed);
        signal = squeeze(acc_healthy(idx, sensor, :));

        ax(1) = nexttile(layout);
        plot(t_healthy, signal, 'k', 'LineWidth', 1);
        title('Sano — 0%');
        ylabel('a [m/s^2]');
        grid on;

        maxAmplitude = max(maxAmplitude, max(abs(signal)));

        for k = 1:numel(severity)
            idx = find_run(meta_dmg, P_fixed, V_fixed, ...
                scenario, severity(k));

            signal = squeeze(acc_dmg(idx, sensor, :));

            ax(k + 1) = nexttile(layout);
            plot(t_dmg, signal, ...
                'Color', colors(k, :), 'LineWidth', 1);

            title(sprintf('Riduzione di rigidezza: %.0f%%', ...
                severity(k) * 100));
            ylabel('a [m/s^2]');
            grid on;

            maxAmplitude = max(maxAmplitude, max(abs(signal)));
        end

        linkaxes(ax, 'xy');
        xlim(ax(1), [0, min(healthy.cfg.t_total, damaged.cfg.t_total)]);

        if maxAmplitude > 0
            ylim(ax(1), 1.05 * [-maxAmplitude, maxAmplitude]);
        end

        xlabel(ax(end), 'Tempo [s]');
        title(layout, sprintf( ...
            ['Effetto del danno — elemento 8 — nodo 10' ...
             ' — %.0f kN/asse — %.0f km/h'], ...
            P_fixed / 1e3, V_fixed));

    otherwise
        error('Caso "%s" non riconosciuto.', index_plot);
end

%% Funzione locale: ricerca della corsa nei metadati
function idx = find_run(meta, loadValue, speedValue, scenario, severity)

    mask = abs(meta(:, 2) - loadValue) < 1e-3 & ...
           abs(meta(:, 3) - speedValue) < 1e-6;

    if nargin == 3
        % Ricerca della condizione sana
        mask = mask & meta(:, 5) == 0;
    else
        % Ricerca della condizione danneggiata
        mask = mask & meta(:, 5) == 1 & ...
               meta(:, 6) == scenario & ...
               abs(meta(:, 7) - severity) < 1e-6;
    end

    idx = find(mask);

    if numel(idx) ~= 1
        error(['Attesa una sola corsa per %.0f N e %.0f km/h; ' ...
               'trovate %d. Controllare dataset, scenario e severità.'], ...
               loadValue, speedValue, numel(idx));
    end
end