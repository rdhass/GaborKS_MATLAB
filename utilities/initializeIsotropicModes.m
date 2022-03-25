function [uhatR,uhatI,vhatR,vhatI,whatR,whatI,kx,ky,kz,...
          gmxloc,gmyloc,gmzloc] ...
          = initializeIsotropicModes(nxQH,nyQH,nzQH,xQH,yQH,zQH,dxQH,dyQH,dzQH,...
          nk,ntheta,kmin,kmax,scalefact,KE,L)
      
    global uhatR uhatI vhatR vhatI whatR whatI
    global kx ky kz
    global gmxloc gmyloc gmzloc

% Step 2: Initialize Gabor modes in each QH region
    % Allocate memory
        gmxloc = zeros(nxQH,nyQH,nzQH,nk*ntheta);
        gmyloc = zeros(nxQH,nyQH,nzQH,nk*ntheta);
        gmzloc = zeros(nxQH,nyQH,nzQH,nk*ntheta);
        kx     = zeros(nxQH,nyQH,nzQH,nk,ntheta);
        ky     = zeros(nxQH,nyQH,nzQH,nk,ntheta);
        kz     = zeros(nxQH,nyQH,nzQH,nk,ntheta);
        uhatR = zeros(nxQH,nyQH,nzQH,nk,ntheta);
        uhatI = zeros(nxQH,nyQH,nzQH,nk,ntheta);
        vhatR = zeros(nxQH,nyQH,nzQH,nk,ntheta);
        vhatI = zeros(nxQH,nyQH,nzQH,nk,ntheta);
        whatR = zeros(nxQH,nyQH,nzQH,nk,ntheta);
        whatI = zeros(nxQH,nyQH,nzQH,nk,ntheta);
    
    % Assign wavevector magnitudes based on logarithmically
    % spaced shells (non-dimensionalized with L)
        kedge = logspace(log10(kmin),log10(kmax),nk+1);
        kmag = (kedge(1:nk)+kedge(2:nk+1))./2;
        dk = kedge(2:nk+1)-kedge(1:nk);

    % Non-dimensional model energy spectrum
        E = getModelSpectrum(kmag,KE,L);
        umag = sqrt(2*E.*dk./ntheta); % Amplitude of each mode such that the sum of the modes gives the correct kinetic energy
    
    for k = 1:nzQH
        zmin = zQH(k);
        for j = 1:nyQH
            ymin = yQH(j);
            for i = 1:nxQH
                xmin = xQH(i);
                
                % Uniformily distribute modes in QH region
                gmxloc(i,j,k,:) = xmin + dxQH*rand(nk*ntheta,1);
                gmyloc(i,j,k,:) = ymin + dyQH*rand(nk*ntheta,1);
                gmzloc(i,j,k,:) = zmin + dzQH*rand(nk*ntheta,1);
                
                for kid = 1:nk
                    % Isotropically sample wave-vector components
                    theta = 2*pi*rand(ntheta,1);
                    kZ = -1 + 2*rand(ntheta,1);
                    
                    % Assign wave-vector componenets
                    kz(i,j,k,kid,:) = kmag(kid)*kZ;
                    r = sqrt(1-kZ.^2);
                    kx(i,j,k,kid,:) = kmag(kid)*r.*cos(theta);
                    ky(i,j,k,kid,:) = kmag(kid)*r.*sin(theta);
                end
            end
        end
    end
    
   
    % Get velocity vector orientations
        % Generate basis of tangent plane to wavevector shell
            p1x = zeros(nxQH,nyQH,nzQH,nk,ntheta);
            p1y = -kz.*ones(nxQH,nyQH,nzQH,nk,ntheta);
            p1z = ky.*ones(nxQH,nyQH,nzQH,nk,ntheta);
            p2x = (ky.^2 + kz.^2).*ones(nxQH,nyQH,nzQH,nk,ntheta);
            p2y = -kx.*ky.*ones(nxQH,nyQH,nzQH,nk,ntheta);
            p2z = -kx.*kz.*ones(nxQH,nyQH,nzQH,nk,ntheta);

            [p1x, p1y, p1z] = normalizeVec(p1x,p1y,p1z);
            [p2x, p2y, p2z] = normalizeVec(p2x,p2y,p2z);
            
        % Assign orientation of velocity vectors
            theta = 2*pi*rand(nxQH,nyQH,nzQH,nk,ntheta);
            orientationX = cos(theta).*p1x + sin(theta).*p2x;
            orientationY = cos(theta).*p1y + sin(theta).*p2y;
            orientationZ = cos(theta).*p1z + sin(theta).*p2z;
            
        % With modulus and orientation in hand we can now assign real and
        % imaginary components
            for kid = 1:nk
                uRmag = rand(1)*sqrt(umag(kid)^2);
                uImag = sqrt(umag(kid)^2 - uRmag.^2);
                
                uhatR(:,:,:,kid,:) = uRmag.*orientationX(:,:,:,kid,:);
                uhatI(:,:,:,kid,:) = uImag.*orientationX(:,:,:,kid,:);
                
                vhatR(:,:,:,kid,:) = uRmag.*orientationY(:,:,:,kid,:);
                vhatI(:,:,:,kid,:) = uImag.*orientationY(:,:,:,kid,:);
                
                whatR(:,:,:,kid,:) = uRmag.*orientationZ(:,:,:,kid,:);
                whatI(:,:,:,kid,:) = uImag.*orientationZ(:,:,:,kid,:);
            end
            
    % Reshape arrays so we just have one big Gabor mode population
    nmodes = nxQH*nyQH*nzQH*nk*ntheta;
    gmxloc = reshape(gmxloc,[nmodes,1]);
    gmyloc = reshape(gmyloc,[nmodes,1]);
    gmzloc = reshape(gmzloc,[nmodes,1]);
    
    kx = reshape(kx,[nmodes,1]);
    ky = reshape(ky,[nmodes,1]);
    kz = reshape(kz,[nmodes,1]);
    
    uhatR = scalefact*reshape(uhatR,[nmodes,1]);
    uhatI = scalefact*reshape(uhatI,[nmodes,1]);
    vhatR = scalefact*reshape(vhatR,[nmodes,1]);
    vhatI = scalefact*reshape(vhatI,[nmodes,1]);
    whatR = scalefact*reshape(whatR,[nmodes,1]);
    whatI = scalefact*reshape(whatI,[nmodes,1]);
    
    % Confirm velocity is divergence free
    [kdotuR,kdotuI] = dotProduct(kx,ky,kz,uhatR,uhatI,vhatR,vhatI,whatR,whatI);
    assert(max(abs(kdotuR)) < 1e-12)
    assert(max(abs(kdotuI)) < 1e-12)
    
% Step3: Strain modes (see appropriate script)
