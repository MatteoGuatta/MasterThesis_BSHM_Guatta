function F = Forze_carichi_viaggianti(t,kd,x_nodes,cro,cri,pos_P_in,V,P_vec)

    F = zeros(kd,1); % forcing vector initialisation (kd: Nr. TOT DOFs)
    pos_P = pos_P_in(:)+V*t; % update on the position of each load
    P_vec = P_vec(:); % vector containing load entities

    idx = find(pos_P >= x_nodes(1) & pos_P <= x_nodes(end)); % Index of active loads on the bridge

    for ip = 1:numel(idx) % cycle over the loads on the bridge

        p = idx(ip); % load index

        xp = pos_P(p); % position of the load

%% Determination of the beam on which the load is

        e = find(x_nodes <= xp,1,'last'); % identification of closest node before the load position
       
        e = min(e,size(cro,1)); % beam element on which the load is now 

%% Determination of the nodes of the beam

        node_1 = cro(e,1); % node 1 of the beam element on which the load is now
       
        node_2 = cro(e,2); % node 2 of the beam element on which the load is now
        
        Le = x_nodes(node_2)-x_nodes(node_1); % length of the beam element on which the load is now
        
        xi = (xp-x_nodes(node_1))/Le; % normalized coordinate of the beam element on which the load is now

%% From loads to nodal forces for the single element

        N = Hermite_EB(xi,Le); % calculation of the shape functions we need for nodal equivalent forces

        fe = -P_vec(p)*N';       % [F1; M1; F2; M2] % nodal equivalent forces

%% Construction of the global forces vector

        for jj = 1:4

            ii = cri(e,jj);
            
            if ii ~= 0 % we need to exclude contrained dofs

                F(ii) = F(ii)+fe(jj);
            
            end
        
        end

    end

 function N = Hermite_EB(xi,Le)
    % I move from concentrated load to nodal equivalent forces and moments
    % w(x,t) = N(xi)*[w1; theta1; w2; theta2].
    N = [1-3*xi^2+2*xi^3, ...
         Le*(xi-2*xi^2+xi^3), ...
         3*xi^2-2*xi^3, ...
         Le*(-xi^2+xi^3)];

 end

end