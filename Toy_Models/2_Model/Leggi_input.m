function fem = Leggi_input(filename)
%READ_BEAM_INP Legge le card di un modello FEM di travi Euler-Bernoulli
% con elementi concentrati
%
%   fem = read_beam_inp('input_file.txt');
%
%   Colonne delle card (unita' SI):
%   *NODES         ID  x  y  vincolo_w  vincolo_theta
%   *BEAMS         ID  nodo1  nodo2  ID_materiale  ID_sezione
%   *MATERIAL      ID  E  rho  alpha  beta
%   *SECTION       ID  A  I
%   *SPRINGS       ID  nodo 1 nodo 2  Kz Ktheta
%   *DAMPERS       ID  nodo 1 nodo 2  cz ctheta
%   *CONC_SPRINGS  ID  nodo  kw  ktheta
%   *CONC_DAMPING  ID  nodo  cw  ctheta
%   *CONC_MASS     ID  nodo  massa  J
%   *END

%   Vincoli: 0 = libero; 1 = bloccato a zero.
%   rho [kg/m^3], A [m^2], I [m^4], J [kg*m^2].
%   kw [N/m], ktheta [N*m/rad], cw [N*s/m], ctheta [N*m*s/rad].
%   alpha [1/s], beta [s], con Ce = alpha*Me + beta*Ke.
%   Le molle e gli smorzatori CONC sono collegati a terra; massa e J
%   sono solidali rispettivamente a w e theta del nodo indicato.

%   OUTPUT: struct contenente le card originali nei campi maiuscoli,
%   piu' i dati pronti per l'assemblaggio:
%     n_nodes, n_el, ndof  numeri di nodi, travi e GDL liberi
%     x_nodes, y_nodes     coordinate, nell'ordine delle righe di NODES
%     bc                   codici originali di vincolo (0/1)
%     cro                  righe di NODES dei due nodi di ogni trave
%     idb                  GDL liberi per nodo: [w theta]; 0 se vincolato
%     cri                  [w1 theta1 w2 theta2] 
%     Le,E,rho,A,I,m       proprieta' per trave; m = rho*A [kg/m]
%     alpha,beta           coefficienti di Rayleigh per trave
%     springs              [kw ktheta] per nodo
%     dampers              [cw ctheta] per nodo
%     nodal_mass           [massa J] per nodo
%
%   Gli ID possono essere non consecutivi e non ordinati. cro contiene
%   INDICI DI RIGA, mentre le card originali mantengono gli ID del file.
%   Le righe CONC sullo stesso nodo vengono sommate. Le card CONC sono
%   facoltative e possono essere vuote. Le altre card devono avere dati.
%   Si accettano righe vuote, commenti % o # e righe che iniziano con **.
%   *END termina la lettura; in sua assenza si legge fino a fine file.

    fid = fopen(filename, 'rt');
    if fid == -1
        error('read_beam_inp:OpenFile', ...
            'Impossibile aprire il file: %s', filename);
    end
    cleanup = onCleanup(@() fclose(fid)); 

    cards = {'NODES', 'BEAMS', 'MATERIAL', 'SECTION', ...
             'CONC_SPRINGS', 'CONC_DAMPING', 'CONC_MASS'....
             'SPRINGS', 'DAMPERS'};

    ncols = [5, 5, 5, 3, 4, 4, 4, 5, 5];

    for k = 1:numel(cards)
        fem.(cards{k}) = zeros(0, ncols(k));
    end

    active = 0;
    line_number = 0;

    %% Lettura: ogni card seleziona la matrice da riempire
    while true
        line = fgetl(fid);
        if ~ischar(line)
            break
        end
        line_number = line_number + 1;
        line = strtrim(regexprep(line, '[%#].*$', ''));
        if isempty(line) || strncmp(line, '**', 2)
            continue
        end

        if line(1) == '*'
            card = upper(strtrim(line(2:end)));
            if strcmp(card, 'END')
                break
            end
            active = find(strcmp(card, cards), 1);
            if isempty(active)
                error('read_beam_inp:UnknownCard', ...
                    'Riga %d: card sconosciuta *%s.', line_number, card);
            end
            continue
        end

        if active == 0
            error('read_beam_inp:MissingCard', ...
                'Riga %d: dati prima di una card.', line_number);
        end

        tokens = regexp(line, '\s+', 'split');
        values = str2double(tokens);
        if numel(values) ~= ncols(active) || ...
                ~isreal(values) || any(~isfinite(values))
            error('read_beam_inp:InvalidRow', ...
                'Riga %d (*%s): attesi %d numeri reali finiti.', ...
                line_number, cards{active}, ncols(active));
        end
        fem.(cards{active})(end+1, :) = values;
    end

    %% Controlli sui dati e sugli identificativi
    for k = 1:4
        if isempty(fem.(cards{k}))
            error('read_beam_inp:MissingData', ...
                'La card *%s e'' assente o vuota.', cards{k});
        end
    end
    for k = 1:numel(cards)
        values = fem.(cards{k});
        ids = values(:, 1);
        if any(ids < 1 | ids ~= fix(ids) | ids > flintmax) || ...
                numel(unique(ids)) ~= numel(ids)
            error('read_beam_inp:InvalidID', ...
                '*%s: gli ID devono essere interi positivi univoci.', cards{k});
        end
    end

    fem.bc = fem.NODES(:, 4:5);
    if any(fem.bc(:) ~= 0 & fem.bc(:) ~= 1)
        error('read_beam_inp:InvalidBC', ...
            '*NODES: i vincoli possono essere solo 0 o 1.');
    end
    if any(any(fem.MATERIAL(:, 2:3) <= 0)) || ...
            any(any(fem.MATERIAL(:, 4:5) < 0))
        error('read_beam_inp:InvalidMaterial', ...
            '*MATERIAL: E e rho devono essere > 0; alpha e beta >= 0.');
    end
    if any(any(fem.SECTION(:, 2:3) <= 0))
        error('read_beam_inp:InvalidSection', ...
            '*SECTION: A e I devono essere > 0.');
    end

    %% Connessioni: conversione degli ID in indici di riga
    fem.n_nodes = size(fem.NODES, 1);
    fem.n_el = size(fem.BEAMS, 1);
    fem.x_nodes = fem.NODES(:, 2);
    fem.y_nodes = fem.NODES(:, 3);
    node_ids = fem.NODES(:, 1);

    fem.cro = lookup_ids(fem.BEAMS(:, 2:3), node_ids, ...
                         '*BEAMS: nodo');
    mat_rows = lookup_ids(fem.BEAMS(:, 4), fem.MATERIAL(:, 1), ...
                          '*BEAMS: materiale');
    sec_rows = lookup_ids(fem.BEAMS(:, 5), fem.SECTION(:, 1), ...
                          '*BEAMS: sezione');

    dx = fem.x_nodes(fem.cro(:, 2)) - fem.x_nodes(fem.cro(:, 1));
    dy = fem.y_nodes(fem.cro(:, 2)) - fem.y_nodes(fem.cro(:, 1));
    fem.Le = hypot(dx, dy);
    if any(fem.Le <= 0) || any(~isfinite(fem.Le))
        error('read_beam_inp:InvalidLength', ...
            '*BEAMS: ogni trave deve avere lunghezza finita e positiva.');
    end

    fem.E     = fem.MATERIAL(mat_rows, 2);
    fem.rho   = fem.MATERIAL(mat_rows, 3);
    fem.alpha = fem.MATERIAL(mat_rows, 4);
    fem.beta  = fem.MATERIAL(mat_rows, 5);
    fem.A     = fem.SECTION(sec_rows, 2);
    fem.I     = fem.SECTION(sec_rows, 3);
    fem.m     = fem.rho .* fem.A;

    %% Numerazione dei soli GDL liberi: w1, theta1, w2, theta2, ...
    fem.idb = zeros(fem.n_nodes, 2);
    fem.ndof = 0;
    for node = 1:fem.n_nodes
        for dof = 1:2
            if fem.bc(node, dof) == 0
                fem.ndof = fem.ndof + 1;
                fem.idb(node, dof) = fem.ndof;
            end
        end
    end
    fem.cri = [fem.idb(fem.cro(:, 1), :), ...
               fem.idb(fem.cro(:, 2), :)];

    %% Elementi concentrati: somma dei contributi per nodo
    fem.ground_springs = collect_nodal(fem.CONC_SPRINGS, node_ids, ...
                                '*CONC_SPRINGS');
    fem.ground_dampers = collect_nodal(fem.CONC_DAMPING, node_ids, ...
                                '*CONC_DAMPING');
     % Collegamenti tra due nodi:
     % colonne = [riga_nodo1, riga_nodo2, valore_w, valore_theta]

     spring_nodes = lookup_ids(fem.SPRINGS(:, 2:3), node_ids, ...
                          '*SPRINGS: nodo');

      damper_nodes = lookup_ids(fem.DAMPERS(:, 2:3), node_ids, ...
                          '*DAMPERS: nodo');

      if any(any(fem.SPRINGS(:, 4:5) < 0))
         error('Leggi_input:NegativeSpring', ...
          '*SPRINGS: le rigidezze devono essere >= 0.');
      end

      if any(any(fem.DAMPERS(:, 4:5) < 0))
         error('Leggi_input:NegativeDamper', ...
          '*DAMPERS: i coefficienti devono essere >= 0.');
      end

      fem.springs = [spring_nodes, fem.SPRINGS(:, 4:5)];
      fem.dampers = [damper_nodes, fem.DAMPERS(:, 4:5)];
   
      fem.nodal_mass = collect_nodal(fem.CONC_MASS, node_ids, ...
                                   '*CONC_MASS');
end

function rows = lookup_ids(ids, available_ids, label)
% Dato un ID del file, restituisce la riga della tabella corrispondente.
    [found, rows] = ismember(ids, available_ids);
    if any(~found(:))
        bad = ids(find(~found, 1));
        error('read_beam_inp:UnknownID', ...
            '%s: ID %g non definito.', label, bad);
    end
end

function values = collect_nodal(card, node_ids, label)
% Una riga per nodo, due colonne: contributo verticale e rotazionale.
    values = zeros(numel(node_ids), 2);
    if isempty(card)
        return
    end
    if any(any(card(:, 3:4) < 0))
        error('read_beam_inp:NegativeValue', ...
            '%s: i valori concentrati devono essere >= 0.', label);
    end
    rows = lookup_ids(card(:, 2), node_ids, [label, ': nodo']);
    for k = 1:size(card, 1)
        values(rows(k), :) = values(rows(k), :) + card(k, 3:4);
    end
end