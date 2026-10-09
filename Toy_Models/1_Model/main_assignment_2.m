%% main:
clear all
close all
clc
%2
%% structure data
%RED BEAM
m_r = 26.2; % [kg/m]
EJ_r = 5.7380e6; % [N*m^2]
EA_r = 6.9076e8; % [N]
L_r = 40; % [m] % 
L_el_r = 4; % [m] 
% GREEN BEAM
m_g = 8.1; % [kg/m]
EJ_g = 3.5400e5; % [N*m^2]
EA_g = 2.1362e8; % [N]
L_g = sqrt(2^2+1.6^2); % [m]
L_el_g = L_g; % [m]
% BLUE BEAM
m_b = 15.8; % [kg/m]
EJ_b = 1.7995e6; % [N*m^2]
EA_b = 4.1586e8; % [N]
L_b = 6; % [m] ?
L_el_b = 2; % [m]
%
fmax = 20; % [Hz]
omegamax = fmax*2*pi; % [rad/s]
eta = 2;% [] safety factor
% ml = ; % [kg] lumped element mass
% Jl = ; % [Kg*m^2] lumped element moment of inertia
% kx = ; % [N/m] spring in x direction
% ky = ; % [N/m] spring in y direction
%% max lenght check
% RED BEAM
L_el_max_r = sqrt((pi^2./(eta*omegamax))*sqrt(EJ_r/m_r));
n_min_r = ceil(L_r/L_el_max_r);
L_el_max_new_r = L_r/n_min_r; % if L_el_max_new > than L_el -> ok
if L_el_max_new_r >= L_el_r
    fprintf('element length of red beam is verified\n')
else
    fprintf('element length of red beam is not verified\n')
end
% GREEN BEAM
L_el_max_g = sqrt((pi^2./(eta*omegamax))*sqrt(EJ_g/m_g));
n_min_g = ceil(L_g/L_el_max_g);
L_el_max_new_g = L_g/n_min_g; % if L_el_max_new > than L_el -> ok
if L_el_max_new_g >= L_el_g
    fprintf('element length of green beam is verified\n')
else
    fprintf('element length of green beam is not verified\n')
end
% BLUE BEAM
L_el_max_b = sqrt((pi^2./(eta*omegamax))*sqrt(EJ_b/m_b));
n_min_b = ceil(L_b/L_el_max_b);
L_el_max_new_b = L_b/n_min_b; % if L_el_max_new > than L_el -> ok
if L_el_max_new_b >= L_el_b
    fprintf('element length of blue beam is verified\n')
else
    fprintf('element length of blue beam is not verified\n')
end


% Loading the *.inp file (file_i) extracting useful information:
%xy: N x 2 matrix containing the coordinates of the nodes
%nnod: number of nodes
%sizee: maximum dimension of the structure
%idb N x 3 matrix numbering each nodal displacement (free and constrained) with different progressive numbers. N is the number of nodes, 3 are the nodal displacements for each node (1, 2, 3) = (x, y, θ)
%ndof: number of free degrees of freedom
%incidenze: N x 6 same idea as for idb, but with N number of elements and 6 the nodal displacement of each element: (1, 2, 3, 4, 5, 6) = (x1, y1, θ 1, x2, y2, θ 2)
%gamma: N x 1 vector containing the orientation of each element with respect to global reference

[file_i,xy,nnod,sizee,idb,ndof,incidenze,l,gamma,m,EA,EJ,posiz,nbeam,pr] = loadstructure;

% Drawing the structure in xy plane
dis_stru(posiz,l,gamma,xy,pr,idb,ndof);


%% Assembling the global mass and stiffness matrices in the global reference system
[M,K] = assem(incidenze,l,m,EA,EJ,gamma,idb);


%% Matrix Partitioning
MFF = M(1:ndof, 1:ndof); % takes only free d.o.f.
KFF = K(1:ndof, 1:ndof);
MCF = M(ndof+1:end, 1:ndof); % takes constrained d.o.f. due to free motion
KCF = K(ndof+1:end, 1:ndof);
MFC = M(1:ndof, ndof+1:end); % takes free d.o.f. due to contrains motion
KFC = K(1:ndof, ndof+1:end);
MCC = M(ndof+1:end, ndof+1:end); % takes contrstrained d.o.f. due to contrains motion
KCC = K(ndof+1:end, ndof+1:end);
%


%% Eigenfrequencies and Mode shapes
[modes, omega2] = eig(KFF,MFF);
omega = sqrt(diag(omega2));
[omega,i_omega] =sort(omega); % puts the omega in order
freq0 = omega/(2*pi);
modes = modes(:,i_omega);


%% plot of modes
n_mode = 2; 
scale_factor = 2;
for i = 1: n_mode
    mode = modes(:,i); % choose the mode you want
    figure(i+1)
    diseg2(mode, scale_factor, incidenze, l, gamma, posiz, idb, xy)
    title(['Mode Shape ', num2str(i), ' - Frequency: ', num2str(freq0(i), '%.2f'), ' Hz'])
    xlabel('[m]')
    ylabel('[m]')
end


%% Damping matrix on the basis of Rayilegh approach
% B = []'; % insert the damping rations of the n modes (identified by modal parameter identification)
% A = zeros(length(B),2);
% for j = 1:length(B)
%     A(j,:) = [1/(2*omega(j)) omega(j)/2];
% end
% ab = A\B; % matrix division to find alpha and beta parameters
ab = [0.1 2e-4;];
C = ab(1)*M + ab(2)*K;

CFF = C(1:ndof, 1:ndof); % top left part of C matrix (free-free)
CCF = C(ndof+1:end, 1:ndof); % bottom left part of C matrix (contrained-free)
CFC = C(1:ndof, ndof+1:end); % top right part of C matrix (free-constrained)
CCC = C(ndof+1:end, ndof+1:end); % bottom right part of C matrix (constrained-constrained)


%% Frequency Response Function (FRF) computation - FEM Approach
% Application of a harmonic force at node 18 in y direction (idb(18,2)) and computation of the FRF at node 10 in y direction (idb(10,2))
F0 = zeros(ndof,1);
F0 (idb(18,2)) = 1000; % The force vector is initialized with zeros and the force is applied at the specified degree of freedom (node 18, y direction).
om = (0:0.01:20)*2*pi; % frequency range (start;discretization;end)
for i = 1:length(om)
    A = -om(i).^2*MFF+1i*om(i)*CFF+KFF;
    X(:,i) = A\F0; % FRF DISPLACEMENT

%    Xpp(:,i) = -om(i)^2*X(:,i); % FRF ACCELLERATION
%    rr(:,i) = A*X(:,i); % FRF FORCE
end
%
frf_index = idb(10,2); % FRF at node 10 in y direction
frf = X(frf_index,:);


%point A
F0 = zeros(ndof,1);
F0 (idb(18,2)) = 1000; % The force vector is initialized with zeros and the force is applied at the specified degree of freedom (node 18, y direction).
om = (0:0.01:20)*2*pi; % frequency range (start;discretization;end)
for i = 1:length(om)
    A = -om(i).^2*MFF+1i*om(i)*CFF+KFF;
    X(:,i) = A\F0; % FRF DISPLACEMENT
%    Xpp(:,i) = -om(i)^2*X(:,i); % FRF ACCELLERATION
%    rr(:,i) = A*X(:,i); % FRF FORCE
end

frf_index = idb(18,2); % FRF at node 18 in y direction
frf = X(frf_index,:);

% figure()
% semilogy(om/(2*pi), abs(X(idb(10,2),:)), 'b', 'LineWidth', 1.5); hold on;
% semilogy(om/(2*pi), abs(X(idb(18,2),:)), 'r', 'LineWidth', 1.5);
% xlabel('f [Hz]')
% ylabel('|FRF| [m/N]')
% title('FRF Amplitude')
% legend('FRF_{B,A}', 'FRF_{A,A}')
% grid on

% figure()
% plot(om/(2*pi), angle(X(idb(10,2),:)), 'b', 'LineWidth', 1.5); hold on;
% plot(om/(2*pi), angle(X(idb(18,2),:)), 'r', 'LineWidth', 1.5);
% xlabel('f [Hz]')
% ylabel('Phase [rad]')
% title('FRF Phase')
% legend('FRF_{B,A}', 'FRF_{A,A}')
% grid on


%
% ndof_R = idb(,)-ndof; % selection of which node and direction for the reaction
% FRF_R = rr(ndof_R ,:);
% figure()
% semilogy(om/(2*pi),abs(FRF_R)); % plot magnitude
% figure()
% plot(om/(2*pi),angle(FRF_R)); % plot phase


%% Frequency Response Function (FRF) computation - Modal Superposition Approach
% Primo modo, Calcolo della matrice di massa, rigidezza e smorzamento modale e della forza modale
Phi_1 = modes(:,1);
Mmod_1 = Phi_1'*MFF*Phi_1;
Kmod_1 = Phi_1'*KFF*Phi_1;
Cmod_1 = Phi_1'*CFF*Phi_1;
Fmod_1 = Phi_1'*F0;

% FRF in modal superposition approach
for ii = 1:length(om)
X_mod_1(:,ii) = (-om(ii)^2*Mmod_1+1i*om(ii)*Cmod_1+Kmod_1)\Fmod_1; % FRF in coordinate modali: 2 righe (modi richiesti) x 2001 colonne (frequenza)
X_m_1(:,ii) = Phi_1*X_mod_1(:,ii); % Passaggio in coordinate reali -> matrice 92 x 2001
end

frf_index_b = idb(10,2); % selection of the node and direction for the frf for point B
frf_index_a = idb(18,2); % selection of the node and direction for the frf for point A

frf_m_1_b = X_m_1(frf_index_b,:);
frf_m_1_a = X_m_1(frf_index_a,:);

% Secondo modo
Phi_2 = modes(:,2);
Mmod_2 = Phi_2'*MFF*Phi_2;
Kmod_2 = Phi_2'*KFF*Phi_2;
Cmod_2 = Phi_2'*CFF*Phi_2;
Fmod_2 = Phi_2'*F0;

%FRF in modal superposition approach
for ii = 1:length(om)
X_mod_2(:,ii) = (-om(ii)^2*Mmod_2+1i*om(ii)*Cmod_2+Kmod_2)\Fmod_2; % FRF in coordinate modali: 2 righe (modi richiesti) x 2001 colonne (frequenza)
X_m_2(:,ii) = Phi_2*X_mod_2(:,ii); % Passaggio in coordinate reali -> matrice 92 x 2001
end
frf_index_b = idb(10,2); % selection of the node and direction for the frf for point B
frf_index_a = idb(18,2); % selection of the node and direction for the frf for point A

frf_m_2_b = X_m_2(frf_index_b,:);
frf_m_2_a = X_m_2(frf_index_a,:);

% Superposition of mode 1  and 2 both for point a and b
frf_m_super_b = frf_m_1_b + frf_m_2_b;
frf_m_super_a = frf_m_1_a + frf_m_2_a;

% semilogy(om/(2*pi),abs(frf_m_1_b)); % plot magnitude
% hold on
% semilogy(om/(2*pi),abs(frf_m_2_b)); % plot magnitude
% figure()
% plot(om/(2*pi),angle(frf_m_1_b)); % plot phase
% hold on
% plot(om/(2*pi),angle(frf_m_2_b)); % plot phase


% Plotting the comparison between the FRF obtained with the FEM approach and the FRF obtained with the modal superposition approach

% figure(10);
% semilogy(om/(2*pi), abs(X(idb(10,2),:)), 'r', 'LineWidth', 1.5);hold on
% semilogy(om/(2*pi),abs(frf_m_super_b), 'ob');
% title('Comparison Fem-Modal FRF Amplitude')
% legend('FRF_{B,A}_{fem}', 'FRF_{B,A}_{modal}')
% ylabel('|FRF| [m/N]')
% xlabel('f [Hz]')
% grid on;
% 
% 
% figure(11);
% semilogy(om/(2*pi), abs(X(idb(18,2),:)), 'r', 'LineWidth', 1.5);hold on
% semilogy(om/(2*pi),abs(frf_m_super_a), 'ob');
% title('Comparison Fem-Modal FRF Amplitude')
% ylabel('|FRF| [m/N]')
% xlabel('f [Hz]')
% legend('FRF_{A,A}_{fem}', 'FRF_{A,A}_{modal}')
% grid on;
% 
% 
% 
% figure(12);
% plot(om/(2*pi), angle(X(idb(10,2),:)), 'r', 'LineWidth', 1.5); hold on;
% semilogy(om/(2*pi),angle(frf_m_super_b), 'ob');
% title('Comparison Fem-Modal FRF Phase')
% xlabel('f [Hz]')
% ylabel('Phase [rad]')
% legend('Phase_{B,A}_{fem}', 'Phase_{B,A}_{modal}')
% grid on;
% 
% figure(13);
% plot(om/(2*pi), angle(X(idb(18,2),:)), 'r', 'LineWidth', 1.5); hold on;
% semilogy(om/(2*pi),angle(frf_m_super_a), 'ob');
% title('Comparison Fem-Modal FRF Phase')
% legend('Phase_{A,A}_{fem}', 'Phase_{A,A}_{modal}')
% xlabel('f [Hz]')
% ylabel('Phase [rad]')
% grid on;


%% Static Response under the weight of the structure

% acceleration vector 
x_dot_dot = zeros(size(xy,1)*3,1);

for i = 2:3:96
    x_dot_dot (i) = -9.81;
end

% computation of gravity force
F_g = M * x_dot_dot;

% Calcolo della forza peso per i soli DoF liberi
F_g_FF = F_g(1:ndof);

% Soluzione del sistema statico solo sui DoF liberi
x_static = KFF \ F_g_FF;


% Inizializza vettore globale degli spostamenti

x_static_full = [x_static; zeros(4,1)];

n_nodi = size(xy,1);        % 32 nodi
disp = zeros(n_nodi, 2);    % [ux, uy] per ogni nodo

for j = 1:n_nodi
    % dof x
    if idb(j,1) > 0
        disp(j,1) = x_static_full(idb(j,1));
    else
        disp(j,1) = 0;
    end

    % dof y
    if idb(j,2) > 0
        disp(j,2) = x_static_full(idb(j,2));
    else
        disp(j,2) = 0;
    end
end

[max_val, max_disp_index] = max(disp(:,2));
fattore_scala = 1000; % per visualizzare meglio
xy_deformato = xy + fattore_scala * (-disp);

% dis_stru2(posiz, l, gamma, xy, pr, idb, ndof, xy_deformato, max_disp_index);


%% Carico Viaggiante 

%% Carico Viaggiante

% 1. Selezione degli elementi appartenenti al Deck (Proprietà 1 - Trave Rossa)
id_prop_deck = 1;
deck_elements_raw = find(pr == id_prop_deck);

% 2. Ordinamento sequenziale lungo la direzione di avanzamento (+X)
% posiz(el, 1) contiene la coordinata X del primo nodo dell'elemento
x_start_el = posiz(deck_elements_raw, 1);
[~, sort_idx] = sort(x_start_el, 'ascend');
deck_elements = deck_elements_raw(sort_idx);

% 3. Parametri cinematici del carico
Mass = 1000;       % [kg]
v_M  = 1;       % [m/s] velocità lungo la via di corsa
dt   = 0.01;   % [s]
dx   = v_M * dt;% [m] avanzamento ad ogni step lungo l'asse dell'elemento
n_modes = 2;
Fn_global = []; % Inizializzazione matrice storie temporali forze globali

% 4. Ciclo di transito sugli elementi del deck ordinati
for idx = 1:length(deck_elements)
    el = deck_elements(idx);
    
    L_el  = l(el);        % Lunghezza dell'elemento corrente (restituita da loadstructure)
    alpha = gamma(el);    % Angolo di inclinazione dell'elemento (restituito da loadstructure)
    
    % Coordinata locale lungo l'asta inclinata
    a = dx:dx:L_el;
    b = L_el - a;
    
    % Scomposizione della forza peso verticale (M*g verso il basso) negli assi locali dell'elemento
    Px_local = -Mass * 9.81 * sin(alpha);
    Py_local = -Mass * 9.81 * cos(alpha);
    
    % Vettore forze nodali locali (Hermite flessione + lineari assiali)
    Fn_local = [
        Px_local * (b ./ L_el);                              % Forza assiale nodo 1
        Py_local * (a.*b.*(b - a)./(L_el^3) + b./L_el);       % Taglio nodo 1
        Py_local * (a.*b.^2 ./ (L_el^2));                    % Momento nodo 1
        Px_local * (a ./ L_el);                              % Forza assiale nodo 2
        Py_local * (a.*b.*(a - b)./(L_el^3) + a./L_el);       % Taglio nodo 2
       -Py_local * (a.^2 .* b ./ (L_el^2))                   % Momento nodo 2
    ];
    
    % Matrice di rotazione da riferimento locale a globale per l'elemento a 2 nodi (6 DoF)
    lam = [ cos(alpha)   sin(alpha)   0;
           -sin(alpha)   cos(alpha)   0;
            0            0            1 ];
    Lambda = blkdiag(lam, lam);
    
    % Rotazione delle forze nodali nel sistema di riferimento globale
    Fn_global_el = Lambda' * Fn_local;
    
    % Assegnazione ai DoF globali dell'elemento
    Temp = zeros(size(M, 1), length(a));
    Temp(incidenze(el, :), :) = Fn_global_el;
    
    % Concatenazione temporale
    Fn_global = [Fn_global Temp];
end

% Aggiunta di 2 secondi di vibrazione libera dopo l'uscita del carico
Fn_global(:, end + round(2/dt)) = 0;

% Vettore temporale completo
t_vect = (1:size(Fn_global, 2)) * dt;

dis_stru(posiz, l, gamma, xy, pr, idb, ndof); hold on;


% Per ciascun elemento della trave del deck (pr == 1):
% Nel file .inp gli elementi del deck collegano i due nodi di estremità.
% Estraiamo le coordinate direttamente da xy:
for k = 1:length(deck_elements)
    el = deck_elements(k);
    
    % Se hai a disposizione la connettività 'incid_nodes' [nodo_1, nodo_2]:
    % plot(xy(incid_nodes(el, 1:2), 1), xy(incid_nodes(el, 1:2), 2), 'y--', 'LineWidth', 3);
    
    % In alternativa, usando posiz (coordinata iniziale), l ed angolo gamma:
    x1 = posiz(el, 1);
    y1 = posiz(el, 2);
    x2 = x1 + l(el) * cos(gamma(el));
    y2 = y1 + l(el) * sin(gamma(el));
    
    plot([x1, x2], [y1, y2], 'y--', 'LineWidth', 3);
end

title('Verifica Deck Selezionato (linea gialla tratteggiata)');


%plotting of the time series of the nodal forces
figure
plot(t_vect, Fn_global);
ylabel('Nodal Forces');
xlabel('time')


% 1. Matrici modali (per i modi selezionati in Phi, ad es. Modi 1 e 2)
Phi = [Phi_1, Phi_2]; % Dimensione: [ndof x n_modes]
n_modes = size(Phi, 2);

% Passaggio a cordinate modali
Qn_global = Phi' * Fn_global(1:ndof, :);

% Plot delle componenti lagrangiane nel tempo
figure;
plot(t_vect, Qn_global(1,:), 'b', 'LineWidth', 1.5); hold on;
plot(t_vect, Qn_global(2,:), 'r--', 'LineWidth', 1.5);
ylabel('Lagrangian Components Q(t) [N]');
xlabel('Time [s]');
legend('Mode 1', 'Mode 2');
grid on;

% Integrazione Numerica con ode45 (Dominio Modale)

Mmod = Phi' * MFF * Phi;
Kmod = Phi' * KFF * Phi;
Cmod = Phi' * CFF * Phi;

% 2. Condizioni iniziali (quiete: velocità e spostamenti nulli)
% Vettore colonna 2*n_modes x 1
X0 = zeros(2 * n_modes, 1);

% 3. Risoluzione del sistema differenziale con ode45
% Notare: t_vect impone a ode45 di restituire i risultati esattamente su quel reticolo temporale
[t_out, X_state] = ode45(@(t, X) odefn(t, X, t_vect, Qn_global, Mmod, Kmod, Cmod), t_vect, X0);

% 4. Ritorno allo spazio fisico
% X_state ha dimensione [Nt x 2*n_modes]
% La prima metà delle colonne sono le velocità modali q_dot, la seconda metà gli spostamenti q
q_displ = X_state(:, n_modes+1:end)'; % Dimensione: [n_modes x Nt]
q_veloc = X_state(:, 1:n_modes)';     % Dimensione: [n_modes x Nt]

% Spostamenti fisici liberi: u(t) = Phi * q(t)
u_phys = Phi * q_displ; % Dimensione: [ndof x Nt]

% 5. Estrazione e plot del segnale su un nodo sensore
% Esempio: spostamento verticale (uy) del nodo 10
nodo_target = 5;
dof_target  = idb(nodo_target, 2); % Grado di libertà verticale

if dof_target > 0 && dof_target <= ndof
    figure;
    plot(t_out, u_phys(dof_target, :) * 1000, 'b', 'LineWidth', 1.5);
    xlabel('Tempo [s]');
    ylabel('Spostamento verticale [mm]');
    title(sprintf('Risposta Dinamica ode45 al Passaggio del Carico - Nodo %d', nodo_target));
    grid on;
end