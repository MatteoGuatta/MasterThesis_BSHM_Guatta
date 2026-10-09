%% main:
clear all
close all
clc
%
%% structure data
m = ; % [kg/m]
EJ = ; % [N*m^2]
EA = ; % [N]
L = ; % [m]
L_el = ; % [m]
fmax = ; % [Hz]
omegamax = fmax*2*pi; % [rad/s]
eta = ;% [] safety factor
ml = ; % [kg] lumped element mass
Jl = ; % [Kg*m^2] lumped element moment of inertia
kx = ; % [N/m] spring in x direction
ky = ; % [N/m] spring in y direction
%% max lenght check
L_el_max = sqrt((pi^2./(eta*omegamax))*sqrt(EJ/m));
n_min = ceil(L/L_el_max);
L_el_max_new = L/n_min; % if L_el_max_new > than L_el -> ok
if L_el_max_new >= L_el
    fprintf('element length is verified')
else
    fprintf('element length is not verified')
end
%% load structure
[file_i,xy,nnod,sizee,idb,ndof,incidenze,l,gamma,m,EA,EJ,posiz,nbeam,pr] = loadstructure;
%% draw structure
dis_stru(posiz,l,gamma,xy,pr,idb,ndof);
%% mass and stiffness matrix
[M,K] = assem(incidenze,l,m,EA,EJ,gamma,idb);
%% matrix partition
MFF = M(1:ndof, 1:ndof); % takes only free d.o.f.
KFF = K(1:ndof, 1:ndof);
MCF = M(ndof+1:end, 1:ndof); % takes constrained d.o.f. due to free motion
KCF = K(ndof+1:end, 1:ndof);
MFC = M(1:ndof, ndof+1:end); % takes free d.o.f. due to contrains motion
KFC = K(1:ndof, ndof+1:end);
MCC = M(ndof+1:end, ndof+1:end); % takes contrstrained d.o.f. due to contrains motion
KCC = K(ndof+1:end, ndof+1:end);
%
%% eigenfrequencies and modeshapes
[modes, omega2] = eig(MFF\KFF);
omega = sqrt(diag(omega2));
[omega,i_omega] =sort(omega); % puts the omega in order
freq0 = omega/(2*pi);
modes = modes(:,i_omega);
%% plot of modes
n_mode = ; % choose how many nodes you want
scale_factor = ;
for i = 1: n_mode
    mode = modes(:,i); % choose the mode you want
    figure()
    diseg2(mode,scale_factor,incidenze,l,gamma,posiz,idb,xy)
end
%% damping matrix
B = []'; % insert the damping rations of the n modes (identified by modal parameter identification)
A = zeros(length(B),2);
for j = 1:length(B)
    A(j,:) = [1/(2*omega(j)) omega(j)/2];
end
ab = A\B; % matrix division to find alpha and beta parameters
C = ab(1)*M+ab(2)*K;
CFF = C(1:ndof, 1:ndof); % top left part of C matrix (free-free)
CCF = K(ndof+1:end, 1:ndof); % bottom left part of C matrix (contrained-free)
CFC = C(1:ndof, ndof+1:end); % top right part of C matrix (free-constrained
CCC = C(ndof+1:end, ndof+1:end); % bottom right part of C matrix (constrained-constrained)
%% frequency responce function
F0 = zeros(ndof,1);
F0 (idb(,)) = 1; % idb is the index node where the force is applied -> idb(node_number, dof_direction)
om = (::)*2*pi; % frequency range (start;discretization;end)
for i = 1:length(om)
    A = -om(i).^2*MFF+1i*om(i)*CFF+KFF;
    X(:,i) = A\F0;
%    Xpp(:,i) = -om(i)^2*X(:,i); % acceleration
%    rr(:,i) = (-om(i).^2*MCF+1i*om(i)*CCF+KCF)*X(:,i); % reaction forces
end
%
frf_index = idb(,); % selection of the node and direction for the frf
frf = X(frf_index,:);
figure()
semilogy(om/(2*pi),abs(frf)); % plot magnitude
figure()
plot(om/(2*pi),angle(frf)); % plot phase
%
ndof_R = idb(,)-ndof; % selection of which node and direction for the reaction
FRF_R = rr(ndof_R ,:);
figure()
semilogy(om/(2*pi),abs(FRF_R)); % plot magnitude
figure()
plot(om/(2*pi),angle(FRF_R)); % plot phase