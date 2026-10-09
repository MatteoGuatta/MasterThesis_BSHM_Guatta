function dis_stru(posiz, l, gamma, xy, pr, idb, ndof)
% DIS_STRU Visualizza la geometria indeformata della struttura, i vincoli e le proprietà delle travi
%
% Input:
%   posiz  - [nbeam x 2] coordinate del primo nodo di ciascun elemento
%   l      - [1 x nbeam] lunghezze degli elementi
%   gamma  - [1 x nbeam] angoli di inclinazione
%   xy     - [nnod x 2] coordinate dei nodi
%   pr     - [1 x nbeam] indici delle proprietà
%   idb    - [nnod x 3] matrice dei gradi di libertà
%   ndof   - numero di GDL liberi

xmax = max(xy(:,1));
xmin = min(xy(:,1));
ymax = max(xy(:,2));
ymin = min(xy(:,2));

dx = (xmax - xmin) / 100;
dy = (ymax - ymin) / 100;
if dx == 0, dx = 0.01; end
if dy == 0, dy = 0.01; end
d = sqrt(dx^2 + dy^2);

% Palette colori per diverse proprietà
colori_palette = [
    1.0  0.0  0.0;   % 1: Rosso (Baseline)
    0.0  0.6  0.2;   % 2: Verde (Danno / Var)
    0.0  0.3  0.9;   % 3: Blu
    0.9  0.5  0.0;   % 4: Arancione
    0.6  0.1  0.7;   % 5: Viola
    0.2  0.2  0.2    % Altro: Grigio scuro
];

colori = zeros(length(pr), 3);
for i = 1:length(pr)
    idx_col = pr(i);
    if idx_col >= 1 && idx_col <= size(colori_palette, 1)
        colori(i, :) = colori_palette(idx_col, :);
    else
        colori(i, :) = [0.2 0.2 0.2];
    end
end

figure('Name', 'Geometria Struttura', 'NumberTitle', 'off');
hold on;

% 1. Disegno degli elementi trave
for i = 1:length(posiz)
    xin = posiz(i, 1);
    yin = posiz(i, 2);
    xfi = posiz(i, 1) + l(i) * cos(gamma(i));
    yfi = posiz(i, 2) + l(i) * sin(gamma(i));
    plot([xin, xfi], [yin, yfi], 'LineWidth', 2, 'Color', colori(i, :));
end
grid on; box on;

% 2. Disegno dei nodi
plot(xy(:, 1), xy(:, 2), 'k.', 'MarkerSize', 15);

% Simboli dei vincoli
triangolo_h = [ 0 0; -sqrt(3)/2  0.5; -sqrt(3)/2 -0.5; 0 0] * d * 2;
triangolo_v = [ 0 0;  0.5 -sqrt(3)/2; -0.5 -sqrt(3)/2; 0 0] * d * 2;
triangolo_r = [ 0 0;  0.5 -sqrt(3)/2; -0.5 -sqrt(3)/2; 0 0] * d * 2 * [sqrt(2)/2 -sqrt(2)/2; -sqrt(2)/2 -sqrt(2)/2];

for ii = 1:size(xy, 1)
    text(xy(ii, 1) + d, xy(ii, 2) + d, num2str(ii), 'FontWeight', 'bold', 'FontSize', 9);
    % Vincolo orizzontale
    if idb(ii, 1) > ndof
        fill(xy(ii, 1) + triangolo_h(:, 1), xy(ii, 2) + triangolo_h(:, 2), 'k');
    end
    % Vincolo verticale
    if idb(ii, 2) > ndof
        fill(xy(ii, 1) + triangolo_v(:, 1), xy(ii, 2) + triangolo_v(:, 2), 'k');
    end
    % Vincolo rotazionale
    if idb(ii, 3) > ndof
        fill(xy(ii, 1) + triangolo_r(:, 1), xy(ii, 2) + triangolo_r(:, 2), 'k');
    end
end

axis equal;
xlabel('X [m]');
ylabel('Y [m]');
title('Geometria Struttura e Condizioni di Vincolo');
end
