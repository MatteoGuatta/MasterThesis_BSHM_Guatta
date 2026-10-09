function modes_plot_conc_element(phi_w, xy_or_fem, freq, moi, incid)
% MODES_PLOT_CONC_ELEMENT Disegna le componenti verticali dei modi di vibrare
% Compatibile sia con struct fem che con matrici standard [xy, incid]
%
% Input:
%   phi_w      : matrice [n_nodi x moi] contenente la componente verticale di ciascun modo
%   xy_or_fem  : matrice coordinate [nnod x 2] o struct fem
%   freq       : vettore delle frequenze proprie [Hz]
%   moi        : numero di modi da plottare
%   incid      : (Opzionale) connettività delle travi [n_el x 2]

if isstruct(xy_or_fem)
    x_nodes = xy_or_fem.x_nodes;
    cro = xy_or_fem.cro;
    n_nodes = xy_or_fem.n_nodes;
    n_el = xy_or_fem.n_el;
else
    x_nodes = xy_or_fem(:, 1);
    n_nodes = size(xy_or_fem, 1);
    if nargin >= 5 && ~isempty(incid)
        cro = incid;
        n_el = size(cro, 1);
    else
        cro = [(1:n_nodes-1)', (2:n_nodes)'];
        n_el = size(cro, 1);
    end
end

beam_nodes = unique(cro(:));
extra_nodes = setdiff((1:n_nodes)', beam_nodes);

for ik = 1:moi
    figure('Name', sprintf('Forma Modale %d', ik), 'NumberTitle', 'off');
    clf;
    hold on;

    % Disegno delle campate / elementi collegati
    for e = 1:n_el
        nodes = cro(e, :);
        plot(x_nodes(nodes), phi_w(nodes, ik), '-o', ...
            'Color', [0 0.4470 0.7410], 'LineWidth', 1.5);
    end

    % Eventuali nodi isolati / TMD / elementi concentrati
    if ~isempty(extra_nodes)
        plot(x_nodes(extra_nodes), phi_w(extra_nodes, ik), 'rs', ...
            'MarkerFaceColor', 'r', 'MarkerSize', 7);
    end

    grid on;
    xlabel('X [m]');
    ylabel(sprintf('\\phi_%d(x)', ik));
    title(sprintf('Modo %d - f = %.3f Hz', ik, freq(ik)));
    hold off;
end
end
