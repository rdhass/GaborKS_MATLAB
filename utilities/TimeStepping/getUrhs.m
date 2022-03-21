function rhs = getUrhs(uhat,dudx,k,Anu,KE,L,nu)
   
    rhs = zeros(3,1);
    delta = eye(3);
    ksq   = sum(k.^2);
    for jj = 1:3
        for kk = 1:3
            for ll = 1:3 % Indices correspond to Pope eqn (11.83)
                rhs(jj) = rhs(jj) - uhat(kk)*dudx(ll,kk)*(delta(jj,ll) - 2*k(jj)*k(ll)/ksq);
            end
        end
    end
    nuk = sqrt(nu^2 + Anu*KE*L*3/4*ksq^(-2/3)) - nu; % Spectral eddy viscosity. Eqn (3.22) in Aditya's thesis
    rhs = rhs - nuk*ksq*uhat;