%% FEM - Lettura input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Integro con Newmark %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DOFs per nodo:
%   1 -> vertical displacement w
%   2 -> rotation theta
%
% Element DOFs:
%   qe = [w1 theta1 w2 theta2]'

clc
clear
close all

%% Folder containing some functions common to other folders

% Cartella che contiene il main in esecuzione
cartella_main = fileparts(mfilename('fullpath'));

% Cartella principale del progetto
cartella_progetto = fileparts(cartella_main);

% Aggiunge la cartella delle funzioni comuni
addpath(fullfile(cartella_progetto, 'Approccio FEM/Funzioni comuni'));

%% Reading of structure data (from the input file)

% dati_modello_fem = Leggi_input('input_file2.txt'); % without TMD
% dati_modello_fem = Leggi_input('input_file3.txt'); % with TMD
dati_modello_fem = Leggi_input('R5_validare_con_Chris.txt'); % Validare con Christian


%% Variables needed for the main to run

n_nodes = dati_modello_fem.n_nodes;
n_el    = dati_modello_fem.n_el;
ndof    = dati_modello_fem.ndof;

x_nodes = dati_modello_fem.x_nodes;
y_nodes = dati_modello_fem.y_nodes;

cro = dati_modello_fem.cro;
idb = dati_modello_fem.idb;
cri = dati_modello_fem.cri;

%% Matrici strutturali - Assembling phase

[M, C, K] = assemble_beam_matrices(dati_modello_fem);

%% Analisi modale 

moi = 7;

[freq, Phi, zeta, modal_par_list] = modal_analysis(M,K,C,moi);

%% Mode shape reconstruction

phi_w = zeros(n_nodes,moi);

for inode = 1:n_nodes

    dof_w = idb(inode,1);

    if dof_w ~= 0

        phi_w(inode,:) = Phi(dof_w,1:moi);

    end

end

%% Plot modi di vibrare

%modes_plot(phi_w,x_nodes,freq,moi)
modes_plot_conc_element(phi_w, dati_modello_fem, freq, moi)

%% Properties of moving loads

V = 43.9;                         % velocita' [m/s]
n_P = 30;                        % numero di carichi
d_P = 5;                       % distanza tra carichi [m]
P_val = 1e5;                     % intensita' di ogni carico [N]

P_vec = P_val*ones(n_P,1);        % anche intensita' diverse, se necessario
pos_P_in = -(0:n_P-1)'*d_P;      % posizioni iniziali [m]

%% Initial conditions for the bridge
kd = ndof;
q0 = zeros(kd,1);                % initial displacements
q_dot0 = zeros(kd,1);            % initial speeds
z0 = [q0; q_dot0];

%% Time interval

t_exit = (x_nodes(end)-min(pos_P_in))/V; % last loads exits
t_free = 1;                      
t_end = t_exit+t_free;

dt_out = 1e-3;                
tspan = linspace(0,t_end,ceil(t_end/dt_out)+1)';

%% Precalcolo matrici

L_M = chol(M,'lower');

MK = -M\K;
  
MC = -M\C;

%% Equivalent nodal forces vector

F_fun = @(tt) Forze_carichi_viaggianti( ...
    tt,kd,x_nodes,cro,cri,pos_P_in,V,P_vec);

%% Parametri di Newmark

gamma_N = 1/2;
beta_N  = 1/4;       % variante delle differenze centrali

%% Passo temporale integrazione

dt = 1e-5;        % passo desiderato [s]
% Griglia uniforme che termina esattamente in t_end
n_steps = ceil(t_end/dt);
t = linspace(0,t_end,n_steps+1)';

%% Integrazione

[q,q_dot,q_ddot] = Newmark_accelerazioni( ...
    M,C,K,F_fun,t,q0,q_dot0,beta_N,gamma_N);

%% Spostamento a mezzeria
% Con n_el = 20, la mezzeria coincide con il nodo 11.

inode_mid = n_el/2 + 1;
dof_mid   = idb(inode_mid,1);

w_mid = q(dof_mid,:)';        % [m], positivo verso l'alto

figure
plot(t,1000*w_mid,'LineWidth',1.5)
grid on

xlabel('Time [s]')
ylabel('Spostamento a mezzeria verso il basso [mm]')
title('Risposta ai carichi viaggianti')
