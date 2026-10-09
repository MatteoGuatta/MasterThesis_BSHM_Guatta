function dis_stru2(posiz,l,gamma,xy,pr,idb,ndof,x_def,max_index)
% Plotta struttura indeformata (blu tratteggiata) e deformata (rossa continua)

xmax = max(xy(:,1)); xmin = min(xy(:,1));
ymax = max(xy(:,2)); ymin = min(xy(:,2));
dx = (xmax - xmin)/100;
dy = (ymax - ymin)/100;
d = sqrt(dx^2 + dy^2);

figure();
hold on;
grid on; box on;
axis equal;

% === STRUTTURA INDEFORMATA (blu tratteggiata) ===
for i=1:length(posiz)
    x0 = posiz(i,:);
    xf = x0 + l(i)*[cos(gamma(i)), sin(gamma(i))];
    plot([x0(1) xf(1)], [x0(2) xf(2)], 'b--', 'LineWidth', 1.5);
end

% === COORDINATE DEFORMATE ===
xy_def = x_def;

% === STRUTTURA DEFORMATA (rossa continua) ===
for i=1:length(posiz)
    x0 = posiz(i,:);
    xf = x0 + l(i)*[cos(gamma(i)), sin(gamma(i))];

    % trova nodi estremi
    node_in = find(all(abs(xy - x0) < 1e-6, 2));
    node_fi = find(all(abs(xy - xf) < 1e-6, 2));
    if isempty(node_in) || isempty(node_fi)
        continue;
    end

    plot([xy_def(node_in,1), xy_def(node_fi,1)], ...
         [xy_def(node_in,2), xy_def(node_fi,2)], 'r-', 'LineWidth', 2);
end

% === PLOT NODI INIZIALI ===
plot(xy_def(:,1), xy_def(:,2), 'ro', 'MarkerFaceColor','r','MarkerSize',6);
plot(xy_def(max_index,1), xy_def(max_index,2) ,'yo', 'MarkerFaceColor','y','MarkerSize',10)
 
% === Vincoli ===
triangolo_h = [ 0 0; -sqrt(3)/2 .5; -sqrt(3)/2 -.5; 0 0]*d*2;
triangolo_v = [ 0 0; .5 -sqrt(3)/2; -.5 -sqrt(3)/2; 0 0]*d*2;
triangolo_r = [0 0; .5 -sqrt(3)/2; -.5 -sqrt(3)/2; 0 0]*d*2 * ...
              [sqrt(2)/2 -sqrt(2/2); -sqrt(2)/2 -sqrt(2)/2];

for ii = 1:size(xy,1)
    text(xy(ii,1) + d, xy(ii,2) + d, num2str(ii));
    if (idb(ii,1) > ndof)
        fill(xy(ii,1) + triangolo_h(:,1), xy(ii,2) + triangolo_h(:,2), 'k');
    end
    if (idb(ii,2) > ndof)
        fill(xy(ii,1) + triangolo_v(:,1), xy(ii,2) + triangolo_v(:,2), 'k');
    end
    if (idb(ii,3) > ndof)
        fill(xy(ii,1) + triangolo_r(:,1), xy(ii,2) + triangolo_r(:,2), 'k');
    end
end


xlabel('[m]')
ylabel('[m]')
% indeformata
h1 = plot([x0(1) xf(1)], [x0(2) xf(2)], 'b--', 'LineWidth', 1.5);

% deformata
h2 = plot([xy_def(node_in,1), xy_def(node_fi,1)], ...
          [xy_def(node_in,2), xy_def(node_fi,2)], 'r-', 'LineWidth', 2);

% nodi deformati
h3 = plot(xy_def(:,1), xy_def(:,2), 'ro', 'MarkerFaceColor','r','MarkerSize',6);

% nodo max spostamento
h4 = plot(xy_def(max_index,1), xy_def(max_index,2), 'yo', ...
          'MarkerFaceColor','y', 'MarkerSize', 10);
legend([h1 h2 h3 h4], ...
       {'Undeformed Structure', ...
        'Deformed Structure', ...
        'Undeformed Nodes', ...
        'Node with maximum displacement'}, ...
        'Location','bestoutside');
