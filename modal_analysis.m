function[freq, Phi, zeta, modal_par_list_table] = modal_analysis(M,K,C,moi)
    
%%%% INPUT %%%%%%%%%%%

    % moi: ordine massimo dei modi di interesse da estrarre
    % M,K,C: matrici strutturali

%%%%%%%%%%%%%%%%%%%%%%

%%%% OUTPUT %%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%


    Ndof_f = size(M);

    if moi > Ndof_f
    
         error('Non hai abbastanza gdl, defiziente!!!')

    end

    modal_par_list = zeros(moi,7);
    
    [Phi,Lambda] = eig(K,M);

    omega2 = diag(Lambda); % squared omega

    omega  = sqrt(omega2); % natural pulsations

    [omega,idx] = sort(omega); % sort the natural pulsations in increasing order

    Phi = Phi(:,idx); % sort the mode shapes according to the increasing order of the natural pulsation

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Normalizzazione modi di vibrare %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for ik = 1:moi

       Phi(:,ik) = Phi(:,ik)./max(abs(Phi(:,ik)));

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    C_mod = Phi'*C*Phi;

    M_mod = Phi'*M*Phi;

    K_mod = Phi'*K*Phi;

    Phi = Phi(:,1:moi);

    freq = omega/(2*pi);

    c_modal_tot = diag(C_mod);

    m_modal_tot = diag(M_mod);

    k_modal_tot = diag(K_mod);
    
    zeta(1:moi) = c_modal_tot(1:moi)./(2.*m_modal_tot(1:moi).*omega(1:moi));

    modal_par_list(1:end,1) = (1:moi)';

    modal_par_list(1:end,2) = m_modal_tot(1:moi);

    modal_par_list(1:end,3) = k_modal_tot(1:moi);

    modal_par_list(1:end,4) = c_modal_tot(1:moi);
   
    modal_par_list(1:end,5) = omega(1:moi);

    modal_par_list(1:end,6) = freq(1:moi);

    modal_par_list(1:end,7) = zeta(1:moi);

    modal_par_list_table = array2table(modal_par_list, ...
    'VariableNames', {'N. ordine','m* [kg]','k* [N/m]','c* [Ns/m]','omega[rad/s]','freq[Hz]','h[-]'});

end