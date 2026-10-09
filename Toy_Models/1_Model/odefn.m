function dX = odefn(t, X, t_vect, Qn_global, Mmod, Kmod, Cmod)
    % Numero di modi considerati
    n_m = size(Mmod, 1);
    
    % X ha dimensione (2*n_m x 1):
    % X(1:n_m)       = velocità modali q_dot
    % X(n_m+1:2*n_m) = spostamenti modali q
    q_dot = X(1:n_m);
    q     = X(n_m+1:end);
    
    % Interpolazione del vettore delle forze modali all'istante continuo t
    Q_t = interp1(t_vect, Qn_global', t, 'linear', 0)'; 
    
    % Calcolo dell'accelerazione modale q_ddot
    q_ddot = Mmod \ (Q_t - Cmod * q_dot - Kmod * q);
    
    % Derivata temporale dello stato: dX = [q_ddot; q_dot]
    dX = [q_ddot; q_dot];
end