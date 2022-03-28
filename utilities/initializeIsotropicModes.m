function [uhatR,uhatI,vhatR,vhatI,whatR,whatI,kx,ky,kz,...
          gmxloc,gmyloc,gmzloc] ...
          = initializeIsotropicModes(nxQH,nyQH,nzQH,xQH,yQH,zQH,dxQH,dyQH,dzQH,...
          nk,ntheta,kmin,kmax,scalefact,KE,L)
  % Function that initializes isotropic modes in each quasi-homogeneous
  % region
  % Inputs:
  %     nxQH, nyQH, nzQH --> Number QH regions in each direction
  %     xQH,yQH,zQH --> QH mesh which includes the boundaries, therefore
  %                     length(xQH) == nxQH+1
  %     dxQH,dyQH,dzQH --> QH mesh grid spacing
  %     nk --> number of wavevector shells per QH region
  %     ntheta --> number of modes per wavevector shell to initialize
  %     kmin,kmax --> The smallest and largest wavenumber modes to
  %                   initialize
  %     scalefact --> Tunable parameter to ensure the induced velocity field has
  %                   the proper energy spectrum 
  %     KE --> Large scale kinetic energy. nxQH X nyQH X nzQH array
  %     L --> integral length scale characteristic of large scale field. 
  %           nxQH X nyQH X nzQH array
      
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
        umag = zeros(nxQH,nyQH,nzQH,nk);
    
    % Assign wavevector magnitudes based on logarithmically
    % spaced shells (non-dimensionalized with L)
        kedge = logspace(log10(kmin),log10(kmax),nk+1);
        kmag = (kedge(1:nk)+kedge(2:nk+1))./2;
        dk = kedge(2:nk+1)-kedge(1:nk);
    
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
                
                % Non-dimensional model energy spectrum
                E = getModelSpectrum(kmag,KE(i,j,k),L(i,j,k));
                
                % Amplitude of each mode such that the sum of the modes gives the correct kinetic energy
                umag(i,j,k,:) = sqrt(2*E.*dk./ntheta);
            end
        end
    end
    
    % Replicate array for fast array operations
    umag = repmat(umag,[1,1,1,1,ntheta]);
   
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
            uRmag = rand(1)*sqrt(umag.^2);
            uImag = sqrt(umag.^2 - uRmag.^2);

            uhatR = uRmag.*orientationX;
            uhatI = uImag.*orientationX;

            vhatR = uRmag.*orientationY;
            vhatI = uImag.*orientationY;

            whatR = uRmag.*orientationZ;
            whatI = uImag.*orientationZ;
            
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
