function modes_plot_conc_element(phi_w, fem, freq, moi)

    % Nodi appartenenti alle travi
    beam_nodes = unique(fem.cro(:));

    % Nodi non appartenenti alle travi, ad esempio il TMD
    extra_nodes = setdiff((1:fem.n_nodes)', beam_nodes);

    for ik = 1:moi

        figure(ik)
        clf
        hold on

        % Disegna soltanto i collegamenti definiti dalle travi
        for e = 1:fem.n_el
            nodes = fem.cro(e, :);

            plot(fem.x_nodes(nodes), phi_w(nodes, ik), ...
                 '-o', 'Color', [0 0.4470 0.7410], ...
                 'LineWidth', 1.5);
        end

        % Nodi aggiuntivi: simboli senza segmenti di collegamento
        if ~isempty(extra_nodes)
            plot(fem.x_nodes(extra_nodes), phi_w(extra_nodes, ik), ...
                 'rs', 'LineStyle', 'none', ...
                 'MarkerFaceColor', 'r', 'MarkerSize', 7);
        end

        grid on
        xlabel('x [m]')
        ylabel(['\phi_', num2str(ik), '(x)'])

        title(['Mode ', num2str(ik), ...
               ' - f = ', num2str(freq(ik), '%.3f'), ' Hz'])

        hold off
    end
end