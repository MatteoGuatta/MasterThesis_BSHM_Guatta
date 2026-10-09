function diseg2(mode, scale_factor, incidenze, l, gamma, posiz, idb, xy)
% DISEG2 Traccia la forma modale deformata della struttura 2D completa
% con interpolazione flessionale ed assiale tramite funzioni di forma.

[n_el, ~] = size(incidenze);
n_gdl = length(mode);

hold on;
for k = 1:n_el
    % Gradi di libertà dell'elemento in coordinate globali
    xkG = zeros(6, 1);
    for iri = 1:6
        if incidenze(k, iri) <= n_gdl && incidenze(k, iri) > 0
            xkG(iri, 1) = mode(incidenze(k, iri));
        else
            xkG(iri, 1) = 0.0;
        end
    end
    
    % Applicazione fattore di scala
    xkG = scale_factor * xkG;
    
    % Rotazione da globale a locale
    lambda = [ cos(gamma(k))  sin(gamma(k)) 0;
              -sin(gamma(k))  cos(gamma(k)) 0;
               0              0             1 ];
    Lambda = blkdiag(lambda, lambda);
    xkL = Lambda * xkG;

    % Calcolo spostamento assiale (u) e trasversale (w) lungo l'elemento
    csi = l(k) * (0:0.05:1);
    
    fu = zeros(6, length(csi));
    fu(1, :) = 1 - csi / l(k);
    fu(4, :) = csi / l(k);
    u = (fu' * xkL)';

    fw = zeros(6, length(csi));
    fw(2, :) = 2 * (csi / l(k)).^3 - 3 * (csi / l(k)).^2 + 1;
    fw(3, :) = l(k) * ((csi / l(k)).^3 - 2 * (csi / l(k)).^2 + csi / l(k));
    fw(5, :) = -2 * (csi / l(k)).^3 + 3 * (csi / l(k)).^2;
    fw(6, :) = l(k) * ((csi / l(k)).^3 - (csi / l(k)).^2);
    w = (fw' * xkL)';

    % Trasformazione da coordinate locali a globali
    xyG   = lambda(1:2, 1:2)' * [u + csi; w];
    undef = lambda(1:2, 1:2)' * [csi; zeros(1, length(csi))];

    % Disegno configurazione indeformata e deformata
    plot(undef(1, :) + posiz(k, 1), undef(2, :) + posiz(k, 2), 'k--', 'LineWidth', 1);
    plot(xyG(1, :)   + posiz(k, 1), xyG(2, :)   + posiz(k, 2), 'b-',  'LineWidth', 1.8);
end

% Nodi deformati
n_nodi = size(idb, 1);
xkG_nodi = zeros(n_nodi, 2);
for k = 1:n_nodi
    for ixy = 1:2
        if idb(k, ixy) <= n_gdl && idb(k, ixy) > 0
            xkG_nodi(k, ixy) = mode(idb(k, ixy));
        end
    end
end
xyG_nodi = xy + scale_factor * xkG_nodi;

plot(xy(:, 1), xy(:, 2), 'k.', 'MarkerSize', 10);
plot(xyG_nodi(:, 1), xyG_nodi(:, 2), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');

grid on;
box on;
axis equal;
end
