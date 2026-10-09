function [M, C, K] = assemble_beam_matrices(fem)
%ASSEMBLE_BEAM_MATRICES Assembla il modello letto da Leggi_input
%
%   fem = read_beam_inp('ponte.inp');
%   [M, C, K] = assemble_beam_matrices(fem);
%
%   Modello: flessione Euler-Bernoulli di travi orizzontali, con
%   w verticale e theta = dw/dx. Massa distribuita consistente.
%   Non comprende deformazione assiale, taglio o inerzia rotatoria
%   distribuita. J in CONC_MASS e' invece un'inerzia nodale aggiuntiva.
%
%   Le matrici sono ndof x ndof, sui soli GDL liberi.
%   I vincoli ideali sono omogenei: i GDL bloccati valgono zero.
%   Le connessioni BEAMS possono essere orientate nei due versi di x.
%
%   Rayleigh viene applicato ALLE SOLE TRAVI, elemento per elemento:
%       Ce = alpha(e)*Me + beta(e)*Ke.
%   Si sommano poi K_conc, C_conc, M_conc. Le molle e masse concentrate
%   non generano automaticamente altri termini di smorzamento Rayleigh.
%   Per un unico materiale: C = alpha*M_beam + beta*K_beam + C_conc.


    %% Controllo travi duplicate tra gli stessi nodi

    %% Controllo di tutti i gruppi di travi duplicate

    if fem.n_el > 1

        % Riconosce come equivalenti anche collegamenti invertiti
        node_pairs = sort(fem.cro, 2);

        [~, ~, group] = unique(node_pairs, 'rows');
        counts = accumarray(group, 1);

        % Tutti i gruppi con almeno due travi sugli stessi nodi
        duplicate_groups = find(counts > 1);

        if ~isempty(duplicate_groups)

            messages = cell(numel(duplicate_groups), 1);

            for ig = 1:numel(duplicate_groups)
                elements = find(group == duplicate_groups(ig));

                % ID originali delle travi
                beam_ids = fem.BEAMS(elements, 1);

                % ID originali dei due nodi comuni
                nodes = node_pairs(elements(1), :);
                node_ids = fem.NODES(nodes, 1);

                messages{ig} = sprintf( ...
                    'Nodi ID %g e ID %g: travi con ID %s.', ...
                    node_ids(1), node_ids(2), mat2str(beam_ids.'));
            end

            error('assemble_beam_matrices:DuplicateBeam', ...
                ['Trovati collegamenti BEAMS duplicati:\n%s\n', ...
                 'Verificare la card *BEAMS.'], ...
                strjoin(messages, sprintf('\n')));
        end
    end

    %% Initialize the matrices

    M = zeros(fem.ndof, fem.ndof);
    C = zeros(fem.ndof, fem.ndof);
    K = zeros(fem.ndof, fem.ndof);

    %% Contributi delle travi
    for e = 1:fem.n_el
        n1 = fem.cro(e, 1);
        n2 = fem.cro(e, 2);
        dx = fem.x_nodes(n2) - fem.x_nodes(n1);
        dy = fem.y_nodes(n2) - fem.y_nodes(n1);
        if abs(dy) > 1e-10*fem.Le(e) || dx == 0
            error('assemble_beam_matrices:InclinedBeam', ...
                ['Trave ID %g: questo assemblatore con GDL [w theta] ', ...
                 'richiede travi orizzontali.'], fem.BEAMS(e, 1));
        end
        Le = fem.Le(e);

        Ke = fem.E(e)*fem.I(e)/Le^3 * ...
            [ 12,     6*Le,  -12,     6*Le;
              6*Le, 4*Le^2, -6*Le, 2*Le^2;
             -12,    -6*Le,   12,    -6*Le;
              6*Le, 2*Le^2, -6*Le, 4*Le^2];

        Me = fem.m(e)*Le/420 * ...
            [156,      22*Le,   54,     -13*Le;
              22*Le,  4*Le^2,  13*Le,  -3*Le^2;
              54,      13*Le,  156,     -22*Le;
             -13*Le, -3*Le^2, -22*Le,   4*Le^2];

        % Se nodo2 e' a sinistra di nodo1, la derivata nella coordinata
        % locale ha segno opposto a theta = dw/dx globale.
        T = diag([1, sign(dx), 1, sign(dx)]);
        Ke = T.' * Ke * T;
        Me = T.' * Me * T;
        Ce = fem.alpha(e)*Me + fem.beta(e)*Ke;

        % cri=0 identifica un GDL bloccato: non compare nel sistema ridotto.
        free = fem.cri(e, :) ~= 0;
        if any(free)
            idx = fem.cri(e, free);
            K(idx, idx) = K(idx, idx) + Ke(free, free);
            M(idx, idx) = M(idx, idx) + Me(free, free);
            C(idx, idx) = C(idx, idx) + Ce(free, free);
        end
    end


    %% Contributi molle e smorzatori a terra e masse concentrate

     for node = 1:fem.n_nodes
        for dof = 1:2
            ii = fem.idb(node, dof);

             if ii ~= 0
                 K(ii, ii) = K(ii, ii) + fem.ground_springs(node, dof);

                 C(ii, ii) = C(ii, ii) + fem.ground_dampers(node, dof);

                 M(ii, ii) = M(ii, ii) + fem.nodal_mass(node, dof);
             end
         end
      end

    %% Contributi tra due nodi

    %% Molle tra due nodi
    for s = 1:size(fem.springs, 1)
        n1 = fem.springs(s, 1);
        n2 = fem.springs(s, 2);

        for dof = 1:2
            % dof = 1: molla verticale
            % dof = 2: molla rotazionale
            k = fem.springs(s, 2 + dof);

            ii = fem.idb(n1, dof);
            jj = fem.idb(n2, dof);

            if ii ~= 0
                K(ii, ii) = K(ii, ii) + k;
            end

            if jj ~= 0
                K(jj, jj) = K(jj, jj) + k;
            end

            if ii ~= 0 && jj ~= 0
                K(ii, jj) = K(ii, jj) - k;
                K(jj, ii) = K(jj, ii) - k;
            end
        end
    end

    %% Smorzatori tra due nodi
    for s = 1:size(fem.dampers, 1)
        n1 = fem.dampers(s, 1);
        n2 = fem.dampers(s, 2);

        for dof = 1:2
            % dof = 1: smorzatore verticale
            % dof = 2: smorzatore rotazionale
            c = fem.dampers(s, 2 + dof);

            ii = fem.idb(n1, dof);
            jj = fem.idb(n2, dof);

            if ii ~= 0
                C(ii, ii) = C(ii, ii) + c;
            end

            if jj ~= 0
                C(jj, jj) = C(jj, jj) + c;
            end

            if ii ~= 0 && jj ~= 0
                C(ii, jj) = C(ii, jj) - c;
                C(jj, ii) = C(jj, ii) - c;
            end
        end
    end

    %% Individuazione dei GDL collegati alle travi

    % Nodi appartenenti ad almeno una trave
    is_beam_node = false(fem.n_nodes, 1);
    is_beam_node(unique(fem.cro(:))) = true;

    % Ogni riga: [nodo1, nodo2, coefficiente_w, coefficiente_theta]
    links = [fem.springs; fem.dampers];

    connected_to_beam = false(fem.n_nodes, 2);

    for il = 1:size(links, 1)
        n1 = links(il, 1);
        n2 = links(il, 2);

        for dof = 1:2

            % Coefficiente nullo: nessun collegamento per questo GDL
            if links(il, 2 + dof) == 0
                continue
            end

            if is_beam_node(n1)
                connected_to_beam(n2, dof) = true;
            end

            if is_beam_node(n2)
                connected_to_beam(n1, dof) = true;
            end
        end
    end

    %% Controlli sui singoli GDL liberi

    dof_names = {'traslazione verticale', 'rotazione'};
    inertia_names = {'massa', 'inerzia rotazionale'};

    for node = 1:fem.n_nodes
        for dof = 1:2
            ii = fem.idb(node, dof);

            % Esclude i GDL vincolati
            if ii == 0
                continue
            end

            no_mass = all(M(ii, :) == 0);
            no_stiffness = all(K(ii, :) == 0);
            no_damping = all(C(ii, :) == 0);

            % Usa l'ID originale del nodo nel messaggio
            node_id = fem.NODES(node, 1);

            if no_mass && no_stiffness && no_damping

                warning('assemble_beam_matrices:InactiveDOF', ...
                    ['Nodo ID %g: il GDL libero di %s non ha ', ...
                     'massa/inerzia, rigidezza ne'' smorzamento. ', ...
                     'Genera una riga e una colonna nulle in M, C e K. ', ...
                     'Verificare collegamenti e vincoli.'], ...
                    node_id, dof_names{dof});

            elseif no_mass && connected_to_beam(node, dof)

                warning('assemble_beam_matrices:MissingMass', ...
                    ['Nodo ID %g: il GDL libero di %s e'' collegato ', ...
                     'a una trave tramite molla e/o smorzatore, ', ...
                     'ma non possiede %s.'], ...
                    node_id, dof_names{dof}, inertia_names{dof});

            end
        end
    end  


        %% Controllo della matrice di rigidezza: warning

    rcK = rcond(full(K));

    if ~isfinite(rcK) || rcK <= eps
        warning('assemble_beam_matrices:SingularK', ...
            ['K e'' singolare o numericamente quasi singolare ', ...
             '(rcond = %.3e). ', ...
             'Verificare eventuali moti rigidi o meccanismi.'], rcK);
    end

    %% Controllo della matrice di massa: arresto se non valida

    [~, pM] = chol(M, 'lower');

    if pM ~= 0
        error('assemble_beam_matrices:InvalidMassMatrix', ...
            ['M non e'' definita positiva. ', ...
             'Impossibile utilizzare l''integratore attuale. ', ...
             'Verificare i GDL senza massa/inerzia ', ...
             'e i GDL liberi inutilizzati.']);
    end

end
