function F = Forze_carichi_viaggianti(t, ndof, xy, incidenze, l, gamma, posiz, deck_elements, pos_P_in, V, P_vec)
% FORZE_CARICHI_VIAGGIANTI Calcola il vettore delle forze nodali equivalenti
% prodotte da un convoglio di carichi concentrati in transito sul deck (2D).
%
% Input:
%   t             : istante temporale corrente [s]
%   ndof          : numero di GDL liberi della struttura
%   xy            : [nnod x 2] coordinate nodali
%   incidenze     : [nbeam x 6] GDL globali degli elementi trave
%   l             : [1 x nbeam] lunghezze elementi
%   gamma         : [1 x nbeam] angoli di inclinazione elementi
%   posiz         : [nbeam x 2] coordinate nodo iniziale elementi
%   deck_elements : vettore indici degli elementi che costituiscono la via di corsa
%   pos_P_in      : vettore posizioni iniziali dei carichi [m]
%   V             : velocità del convoglio [m/s]
%   P_vec         : intensità dei carichi [N] (positivi verso il basso)
%
% Output:
%   F             : [ndof x 1] vettore forze nodali equivalenti sui GDL liberi

F = zeros(ndof, 1);

pos_P = pos_P_in(:) + V * t;
P_vec = P_vec(:);

% Coordinate X di inizio e fine del deck
deck_el      = deck_elements(:);
x_start_deck = posiz(deck_el, 1);
l_deck       = l(deck_el); l_deck = l_deck(:);
gamma_deck   = gamma(deck_el); gamma_deck = gamma_deck(:);
x_end_deck   = x_start_deck + l_deck .* cos(gamma_deck);

x_deck_start = min([x_start_deck; x_end_deck]);
x_deck_end   = max([x_start_deck; x_end_deck]);

% Identifica i carichi attualmente presenti sul ponte
active_p = find(pos_P >= x_deck_start - 1e-9 & pos_P <= x_deck_end + 1e-9);

for ip = 1:numel(active_p)
    p  = active_p(ip);
    xp = pos_P(p);
    P_val = P_vec(p);

    % Ricerca dell'elemento del deck su cui si trova il carico xp
    el_found = 0;
    for idx = 1:numel(deck_elements)
        el = deck_elements(idx);
        x1 = posiz(el, 1);
        x2 = x1 + l(el) * cos(gamma(el));
        
        x_min = min(x1, x2);
        x_max = max(x1, x2);
        
        if xp >= x_min && xp <= x_max
            el_found = el;
            break;
        end
    end

    if el_found == 0
        continue;
    end

    el = el_found;
    L_el  = l(el);
    alpha = gamma(el);

    % Coordinata locale lungo l'asta
    dx_local = xp - posiz(el, 1);
    if abs(cos(alpha)) > 1e-6
        a = dx_local / cos(alpha);
    else
        a = 0;
    end
    a = max(0, min(L_el, a));
    b = L_el - a;
    xi = a / L_el;

    % Scomposizione del carico gravitazionale (diretto verso il basso -Y globale)
    Px_local = -P_val * sin(alpha);
    Py_local = -P_val * cos(alpha);

    % Funzioni di forma Hermite (flessione) e lineari (assiale)
    % GDL locali: [u1; w1; theta1; u2; w2; theta2]
    fe_local = [
        Px_local * (b / L_el);                                    % Assiale nodo 1
        Py_local * (1 - 3*xi^2 + 2*xi^3);                         % Taglio nodo 1
        Py_local * (L_el * (xi - 2*xi^2 + xi^3));                 % Momento nodo 1
        Px_local * (a / L_el);                                    % Assiale nodo 2
        Py_local * (3*xi^2 - 2*xi^3);                             % Taglio nodo 2
        Py_local * (L_el * (-xi^2 + xi^3))                        % Momento nodo 2
    ];

    % Matrice di rotazione locale -> globale
    lam = [ cos(alpha),  sin(alpha), 0;
           -sin(alpha),  cos(alpha), 0;
            0,           0,          1 ];
    Lambda = blkdiag(lam, lam);

    % Forze nodali equivalenti nel sistema globale
    fe_global = Lambda' * fe_local;

    % Assegnazione ai gradi di libertà liberi della struttura
    for jj = 1:6
        ii = incidenze(el, jj);
        if ii > 0 && ii <= ndof
            F(ii) = F(ii) + fe_global(jj);
        end
    end
end

end
