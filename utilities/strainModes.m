function strainModes(xLES,yLES,zLES,xF,yF,zF,gradU,U,V,W,gmxloc,gmyloc,gmzloc,...
    kx,ky,kz,uhatR,uhatI,vhatR,vhatI,whatR,whatI,Anu,KE,L,nu,ctau)
    
    % Tell MATLAB to modify global variables
        global uhatR vhatR whatR;
        global uhatI vhatI whatI;
        global kx ky kz;
        global gmxloc gmyloc gmzloc;
        
    % Compute eddy lifetime
        nmodes = length(gmxloc);
        S = sqrt(squeeze(sum(sum(gradU.*gradU,1),2)));
        Sm = interp3(xLES,yLES,zLES,S,gmxloc,gmyloc,gmzloc,'spline');
        kabs = sqrt(kx.^2 + ky.^2 + kz.^2);    
        tau = ctau./Sm.*kabs.^(-2/3)./sqrt(hypergeom([1/3,17/6],4/3,-1./kabs.^2));
    %     save('./data/EddyLifetime.mat','tau','S','kabs')


    % Interpolate gradU to mode location
        dudxGM = zeros(3,3,nmodes);
        for j = 1:3
            for i = 1:3
                dudxGM(i,j,:) = interp3(xLES,yLES,zLES,squeeze(gradU(i,j,:,:,:)),gmxloc,gmyloc,gmzloc,'spline');
            end
        end

    % Determine minimum stable time step
        umax = max(sqrt(uhatR.^2 + uhatI.^2 + vhatR.^2 + vhatI.^2 + whatR.^2 + whatI.^2));
        Umax = max(sqrt(U(:).^2 + V(:).^2 + W(:).^2));
        dxF = xF(2)-xF(1);
        dt = min([dxF/umax,dxF/Umax,1/max(S(:))])/4/2;
        
        pool = parpool('threads');
        parfor n = 1:nmodes
            k = [kx(n); ky(n); kz(n)];
            uR = [uhatR(n); vhatR(n); whatR(n)];
            uI = [uhatI(n); vhatI(n); whatI(n)];
            dudx = dudxGM(:,:,n);

            for tid = 1:round(tau(n)/dt)
                [uR, uI, k] = rk4Step(uR,uI,k,dt,Anu,KE,L,nu,dudx);
    %             assert(max(abs(k)) < 1e2);
                assert(max(abs(uR)) < 1e2);
                assert(max(abs(uI)) < 1e2);
            end

            kx(n) = k(1); ky(n) = k(2); kz(n) = k(3);
            uhatR(n) = uR(1); vhatR(n) = uR(2); whatR(n) = uR(3);
            uhatI(n) = uI(1); vhatI(n) = uI(2); whatI(n) = uI(3); 
            if mod(n,500) == 0; disp([num2str(n/nmodes*100),'% complete']); end
        end
        delete(pool)
