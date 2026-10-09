function [q,q_dot,q_ddot] = Newmark_accelerazioni( ...
    M,C,K,F_fun,t,q0,q_dot0,beta_N,gamma_N)

    % The integrating step dt is constant (different from the download dt,
    % called dt_out)

    dt = t(2)-t(1);

    kd = size(M,1);
    nt = numel(t);

    %% Initilize vectors

    q      = zeros(kd,nt); % position
    q_dot  = zeros(kd,nt); % velocity
    q_ddot = zeros(kd,nt); % acceleration

    q(:,1)     = q0(:);
    q_dot(:,1) = q_dot0(:);

    % Accelerazione iniziale dall'equilibrio dinamico
    q_ddot(:,1) = M\(F_fun(t(1)) - C*q_dot(:,1) - K*q(:,1));

    %% Matrice efficace per le accelerazioni

    A_eff = M + gamma_N*dt*C + beta_N*dt^2*K;

    L_eff = chol(A_eff,'lower');

    %% Integrazione temporale

    prossima_percentuale = 1;
    fprintf('Simulation status: 0%%\n');

    for n = 1:nt-1

        % Predittori
        q_pred = q(:,n) + dt*q_dot(:,n) ...
            + dt^2*(0.5-beta_N)*q_ddot(:,n);

        v_pred = q_dot(:,n) ...
            + dt*(1-gamma_N)*q_ddot(:,n);

        % Forze al passo successivo
        F_next = F_fun(t(n+1));

        % Accelerazione al passo successivo
        rhs = F_next - C*v_pred - K*q_pred;

        q_ddot(:,n+1) = L_eff'\(L_eff\rhs);

        % Aggiornamento di spostamenti e velocita'
        q(:,n+1) = q_pred + beta_N*dt^2*q_ddot(:,n+1);

        q_dot(:,n+1) = v_pred + gamma_N*dt*q_ddot(:,n+1);


        % Percentuale di integrazione completata
        percentuale = floor(100*n/(nt-1));

        % Notifica ogni soglia dell'1% raggiunta
        while prossima_percentuale <= percentuale

           fprintf('Simulation status: %d%%\n', prossima_percentuale);

           prossima_percentuale = prossima_percentuale + 1;

        end

end