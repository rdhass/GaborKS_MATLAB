function strainModes(xLES,yLES,zLES,xQHcent,yQHcent,zQHcent,gradU,U,V,W,gmxloc,gmyloc,gmzloc,...
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
    
    % Specify L and KE for the mode location. This is required for the
    % spectral eddy viscosity model
    Lgm = interp3(xQHcent',yQHcent,zQHcent,L,gmxloc,gmyloc,gmzloc,'spline');
    KEgm = interp3(xQHcent',yQHcent,zQHcent,KE,gmxloc,gmyloc,gmzloc,'spline');

    % Find max velocities for determining minimum stable time step below
    umax = max(sqrt(uhatR.^2 + uhatI.^2 + vhatR.^2 + vhatI.^2 + whatR.^2 + whatI.^2));
    Umax = max(sqrt(U(:).^2 + V(:).^2 + W(:).^2));
%     dt = min([min(1./(kabs*umax)),min(1./(kabs*Umax)),1/max(S(:))])/4/2;
        
        pool = parpool('threads');
        parfor n = 1:nmodes
            k = [kx(n); ky(n); kz(n)];
            uR = [uhatR(n); vhatR(n); whatR(n)];
            uI = [uhatI(n); vhatI(n); whatI(n)];
            dudx = dudxGM(:,:,n);
            kmag = sqrt(sum(k.^2));
            dt = min(1/(kmag*umax),1/max(S(:)))/4/2;
            
            for tid = 1:round(tau(n)/dt)
                [uR, uI, k] = rk4Step(uR,uI,k,dt,Anu,KEgm(n),Lgm(n),nu,dudx);
                assert(max(abs(uR)) < 1e2);
                assert(max(abs(uI)) < 1e2);
            end

            kx(n) = k(1); ky(n) = k(2); kz(n) = k(3);
            uhatR(n) = uR(1); vhatR(n) = uR(2); whatR(n) = uR(3);
            uhatI(n) = uI(1); vhatI(n) = uI(2); whatI(n) = uI(3); 
            if mod(n,500) == 0; disp([num2str(n/nmodes*100),'% complete']); end
        end
        delete(pool)
