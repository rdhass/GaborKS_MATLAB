function [uRout, uIout, kout] = rk4StepV1(uR,uI,k,dt,Anu,epsilon,dudx)
    % RK4 coefficients
        alpha = [0 0.5 0.5 1];
        beta  = [1/6 1/3 1/3 1/6];
    
    delta = eye(3);
    duR = zeros(3,1);
    duI = zeros(3,1);
    dk = zeros(3,1);
    rhsR = zeros(3,1);
    rhsI = zeros(3,1);
    rhsK = zeros(3,1);
    for rk = 1:4
        uRstar = uR + dt*alpha(rk)*rhsR;
        uIstar = uI + dt*alpha(rk)*rhsI;
        kstar  = k  + dt*alpha(rk)*rhsK;
        ksq   = sum(kstar.^2);
        for jj = 1:3
            for kk = 1:3
                for ll = 1:3 % Indices correspond to Pope eqn (11.83)
                    rhsR(jj) = rhsR(jj) - uRstar(kk)*dudx(ll,kk)*(delta(jj,ll) - 2*kstar(jj)*kstar(ll)/ksq^2);
                    rhsI(jj) = rhsI(jj) - uIstar(kk)*dudx(ll,kk)*(delta(jj,ll) - 2*kstar(jj)*kstar(ll)/ksq^2);
                end
                rhsK(jj) = rhsK(jj) - kstar(kk)*dudx(kk,jj);
            end
        end
        nuk = Anu*epsilon^(1/3)*ksq^(-2/3);
        rhsR = rhsR - nuk*ksq*uRstar;
        rhsI = rhsI - nuk*ksq*uIstar;

        duR = duR + dt*beta(rk)*rhsR;
        duI = duI + dt*beta(rk)*rhsI;
        dk  = dk  + dt*beta(rk)*rhsK;
    end
    uRout = uR + duR;
    uIout = uI + duI;
    kout  = k  + dk;