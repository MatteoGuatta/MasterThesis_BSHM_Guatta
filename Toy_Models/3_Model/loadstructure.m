function [file_i, xy, nnod, sizee, idb, ndof, incidenze, l, gamma, m, EA, EJ, posiz, nbeam, pr, conc_data] = loadstructure(file_i)
% LOADSTRUCTURE Carica la geometria e le proprietà della struttura da file *.inp
% Supporta sia chiamate interattive che automatiche (passando il nome del file).
%
% Output:
%   file_i     - Nome del file caricato
%   xy         - [nnod x 2] Coordinate nodali (x, y)
%   nnod       - Numero di nodi
%   sizee      - Dimensione massima caratteristica della struttura
%   idb        - [nnod x 3] Matrice dei gradi di libertà (prima liberi 1..ndof, poi vincolati)
%   ndof       - Numero di GDL liberi
%   incidenze  - [nbeam x 6] GDL globali associati a ciascun elemento
%   l          - [1 x nbeam] Lunghezza di ciascun elemento
%   gamma      - [1 x nbeam] Angolo di orientazione di ciascun elemento (rad)
%   m          - [1 x nbeam] Massa per unità di lunghezza [kg/m]
%   EA         - [1 x nbeam] Rigidezza assiale [N]
%   EJ         - [1 x nbeam] Rigidezza flessionale [N*m^2]
%   posiz      - [nbeam x 2] Coordinate del primo nodo di ciascun elemento
%   nbeam      - Numero totale di elementi trave
%   pr         - [1 x nbeam] Indice di proprietà assegnato a ciascuna trave
%   conc_data  - Struct con elementi concentrati (molle, masse, smorzatori, deck)

if nargin < 1 || isempty(file_i)
    disp(' ');
    file_i = input(' Name of the input file *.inp (without extension) = ', 's');
    disp(' ');
end

% Gestione estensione .inp
[fpath, fname, fext] = fileparts(file_i);
if isempty(fext)
    file_full = fullfile(fpath, [fname, '.inp']);
else
    file_full = file_i;
    file_i = fullfile(fpath, fname);
end

if exist(file_full, 'file') ~= 2
    error(['Il file di input non esiste: ', file_full]);
end

fid_i = fopen(file_full, 'r');
if fid_i == -1
    error(['Impossibile aprire il file: ', file_full]);
end

%% 1. Lettura card *NODES
findcard(fid_i, '*NODES', false);
iconta = 0;
while true
    line = scom(fid_i);
    if isempty(line) || ~isempty(strfind(line, '*ENDNODES'))
        break;
    end
    tmp = sscanf(line, '%i %i %i %i %f %f')';
    if length(tmp) >= 6
        iconta = iconta + 1;
        if iconta ~= tmp(1)
            error('Errore: nodi non numerati in ordine progressivo');
        end
        ivinc(iconta, :) = tmp(2:4);
        xy(iconta, :)    = tmp(5:6);
    end
end
nnod = iconta;
sizee = sqrt((max(xy(:,1)) - min(xy(:,1)))^2 + (max(xy(:,2)) - min(xy(:,2)))^2);

%% 2. Costruzione matrice IDB (prima GDL liberi 1..ndof, poi vincolati)
idof = 0;
idb = zeros(nnod, 3);
for i = 1:nnod
    for j = 1:3
        if ivinc(i, j) == 0
            idof = idof + 1;
            idb(i, j) = idof;
        end
    end
end
ndof = idof;

for i = 1:nnod
    for j = 1:3
        if ivinc(i, j) == 1
            idof = idof + 1;
            idb(i, j) = idof;
        end
    end
end

%% 3. Lettura card *PROPERTIES
findcard(fid_i, '*PROPERTIES', false);
iconta = 0;
properties = [];
while true
    line = scom(fid_i);
    if isempty(line) || ~isempty(strfind(line, '*ENDPROPERTIES'))
        break;
    end
    tmp = sscanf(line, '%i %f %f %f')';
    if length(tmp) >= 4
        iconta = iconta + 1;
        properties(tmp(1), :) = tmp(2:4); % [m, EA, EJ]
    end
end

%% 4. Lettura card *BEAMS
findcard(fid_i, '*BEAMS', false);
iconta = 0;
incid = [];
while true
    line = scom(fid_i);
    if isempty(line) || ~isempty(strfind(line, '*ENDBEAMS'))
        break;
    end
    tmp = sscanf(line, '%i %i %i %i')';
    if length(tmp) >= 4
        iconta = iconta + 1;
        incid(iconta, :) = tmp(2:3);
        pr(1, iconta)    = tmp(4);
        m(1, iconta)     = properties(tmp(4), 1);
        EA(1, iconta)    = properties(tmp(4), 2);
        EJ(1, iconta)    = properties(tmp(4), 3);
        
        n1 = incid(iconta, 1);
        n2 = incid(iconta, 2);
        dx = xy(n2, 1) - xy(n1, 1);
        dy = xy(n2, 2) - xy(n1, 2);
        
        l(1, iconta)     = sqrt(dx^2 + dy^2);
        gamma(1, iconta) = atan2(dy, dx);
        incidenze(iconta, :) = [idb(n1, :), idb(n2, :)];
        posiz(iconta, :) = xy(n1, :);
    end
end
nbeam = iconta;

%% 5. Card Opzionali per Elementi Concentrati e Deck
conc_data.springs = [];  % [nodo, kx, ky, ktheta]
conc_data.masses  = [];  % [nodo, mx, my, Jtheta]
conc_data.dampers = [];  % [nodo, cx, cy, ctheta]
conc_data.deck    = [];  % [id_beam1, id_beam2, ...]

% A. *CONC_SPRINGS
if findcard(fid_i, '*CONC_SPRINGS', true)
    while true
        line = scom(fid_i);
        if isempty(line) || ~isempty(strfind(line, '*ENDCONC_SPRINGS'))
            break;
        end
        tmp = sscanf(line, '%i %f %f %f')';
        if length(tmp) >= 4
            conc_data.springs = [conc_data.springs; tmp];
        end
    end
end

% B. *CONC_MASS
if findcard(fid_i, '*CONC_MASS', true)
    while true
        line = scom(fid_i);
        if isempty(line) || ~isempty(strfind(line, '*ENDCONC_MASS'))
            break;
        end
        tmp = sscanf(line, '%i %f %f %f')';
        if length(tmp) >= 4
            conc_data.masses = [conc_data.masses; tmp];
        end
    end
end

% C. *CONC_DAMPERS
if findcard(fid_i, '*CONC_DAMPERS', true)
    while true
        line = scom(fid_i);
        if isempty(line) || ~isempty(strfind(line, '*ENDCONC_DAMPERS'))
            break;
        end
        tmp = sscanf(line, '%i %f %f %f')';
        if length(tmp) >= 4
            conc_data.dampers = [conc_data.dampers; tmp];
        end
    end
end

% D. *DECK (elementi su cui transitano i carichi)
if findcard(fid_i, '*DECK', true)
    while true
        line = scom(fid_i);
        if isempty(line) || ~isempty(strfind(line, '*ENDDECK'))
            break;
        end
        tmp = sscanf(line, '%i')';
        if ~isempty(tmp)
            conc_data.deck = [conc_data.deck, tmp];
        end
    end
end

fclose(fid_i);
end
