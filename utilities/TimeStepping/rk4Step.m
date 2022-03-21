function [uRout, uIout, kout] = rk4Step(uR,uI,k,dt,Anu,KE,L,nu,dudx)
    % RK4 coefficients
        alpha = [0 0.5 0.5 1];
        beta  = [1/6 1/3 1/3 1/6];
    
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
        rhsR = getUrhs(uRstar,dudx,kstar,Anu,KE,L,nu);
        rhsI = getUrhs(uIstar,dudx,kstar,Anu,KE,L,nu);
        rhsK = getKrhs(kstar,dudx);

        duR = duR + dt*beta(rk)*rhsR;
        duI = duI + dt*beta(rk)*rhsI;
        dk  = dk  + dt*beta(rk)*rhsK;
    end
    uRout = uR + duR;
    uIout = uI + duI;
    kout  = k  + dk;
