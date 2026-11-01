function found = findcard(fid, card, is_optional)
% FINDCARD Cerca una card nel file specificato da fid.
% Posiziona il puntatore sulla riga successiva alla card.
% Restituisce 1 se trovata, 0 se opzionale e non trovata.

if nargin < 3
    is_optional = false;
end

maxiter = 1e5;
frewind(fid);

card_clean = upper(strtrim(card));

for i = 1:maxiter
    if feof(fid)
        if is_optional
            found = 0;
            frewind(fid);
            return;
        else
            error(['Card obbligatoria non trovata nel file: ', card]);
        end
    end

    riga = fgets(fid);
    riga_trim = upper(strtrim(riga));

    % Ignora righe vuote o commenti
    if isempty(riga_trim) || riga_trim(1) == '!'
        continue;
    end

    % Se stiamo cercando una card che inizia con '*', non dobbiamo matchare '*END...'
    if startsWith(card_clean, '*') && ~startsWith(card_clean, '*END')
        if startsWith(riga_trim, '*END')
            continue;
        end
    end

    % Verifica corrispondenza con la card cercata
    tokens = strsplit(riga_trim);
    if ~isempty(tokens) && strcmp(tokens{1}, card_clean)
        found = 1;
        return;
    end
end

if is_optional
    found = 0;
    frewind(fid);
else
    error(['Card obbligatoria non trovata nel file: ', card]);
end
