function [uout,vout,wout] = renderVelocityXYperiodic(uhatR,uhatI,vhatR,vhatI,whatR,whatI,...
    kx,ky,kz,gmxloc,gmyloc,gmzloc,nxsupp,nysupp,nzsupp,xFp,yFp,zF)
    % Inputs:
    %   uhatR, uhatI, ... --> complex velocity amplitudes associated with
    %                         each Gabor mode. Array size: nmodesX1
    %   kx, ky, kz --> wave-vector components for each GM. Array size: nmodesX1
    %   gmxloc, gmyloc, gmzloc --> Physical location of each GM. Array size: nmodesX1
    %   nxsupp, nysupp, nzsupp --> Number of fine grid points (on each side of a given GM) 
    %                              spanned by the window function. Scalar
    %                              (integer)
    %   xF, yF, zF --> Eulerian grid for fine mesh. Array size: 1XnF
    % Outputs:
    %   u, v, w --> Physical space velocities. Array size: nxF X nyF X nzF
        
    global uhatR uhatI vhatR vhatI whatR whatI
    global kx ky kz
    global gmxloc gmyloc gmzloc
    
    dxF = xFp(2)-xFp(1);
    dyF = yFp(2)-yFp(1);
    dzF = zF(2)-zF(1);
    
    nxFp = length(xFp);
    nyFp = length(yFp);
    
    nxF = nxFp - nxsupp;
    nyF = nyFp - nysupp;
    nzF = length(zF);
        
    u = zeros(nxFp,nyFp,nzF);
    v = u;
    w = u;
    
    [xFp,yFp,zF] = ndgrid(xFp,yFp,zF);
    
    wxSupport = (nxsupp + 1)*dxF;
    wySupport = (nysupp + 1)*dyF;
    wzSupport = (nzsupp + 1)*dzF;
    
    nmodes = length(kx);
 
    % Given a mode location and its support width we can convert this
    % information to indices on the fine grid 
    for n = 1:nmodes
        ist = ceil(gmxloc(n)/dxF);%  - nxsupp/2;
        ien = floor(gmxloc(n)/dxF) + nxsupp;%/2;
        
        jst = ceil(gmyloc(n)/dyF);%  - nysupp/2;
        jen = floor(gmyloc(n)/dyF) + nysupp;%/2;
        
        kst = max([1,ceil(gmzloc(n)/dzF)  - nzsupp/2]);
        ken = min([nzF,floor(gmzloc(n)/dzF) + nzsupp/2]);
        
        x = xFp(ist:ien,jst:jen,kst:ken);
        y = yFp(ist:ien,jst:jen,kst:ken);
        z = zF(ist:ien,jst:jen,kst:ken);
        
        kdotx = kx(n).*(x-gmxloc(n)) + ky(n).*(y-gmyloc(n)) + kz(n).*(z-gmzloc(n));
        cs = cos(kdotx);
        ss = sin(kdotx);
        
        % Window function
        fx = cos(pi*(x - gmxloc(n))/wxSupport);
        fy = cos(pi*(y - gmyloc(n))/wySupport);
        fz = cos(pi*(z - gmzloc(n))/wzSupport);
        f = fx.*fy.*fz;
        
        % Resulting velocity field
        u(ist:ien,jst:jen,kst:ken) = u(ist:ien,jst:jen,kst:ken) + f.*(2*uhatR(n).*cs - 2*uhatI(n).*ss);
        v(ist:ien,jst:jen,kst:ken) = v(ist:ien,jst:jen,kst:ken) + f.*(2*vhatR(n).*cs - 2*vhatI(n).*ss);
        w(ist:ien,jst:jen,kst:ken) = w(ist:ien,jst:jen,kst:ken) + f.*(2*whatR(n).*cs - 2*whatI(n).*ss);
        
        if mod(n,1000) == 0
          disp([num2str(n/nmodes*100),'% Complete'])
        end
    end
    
     % Add periodic contribution
    u(nxsupp/2+1:nxsupp,:,:) = u(nxsupp/2+1:nxsupp,:,:) + u(nxFp-nxsupp/2+1:nxFp,:,:);
    u(:,nysupp/2+1:nysupp,:) = u(:,nysupp/2+1:nysupp,:) + u(:,nyFp-nysupp/2+1:nyFp,:);
    
    v(nxsupp/2+1:nxsupp,:,:) = v(nxsupp/2+1:nxsupp,:,:) + v(nxFp-nxsupp/2+1:nxFp,:,:);
    v(:,nysupp/2+1:nysupp,:) = v(:,nysupp/2+1:nysupp,:) + v(:,nyFp-nysupp/2+1:nyFp,:);
    
    w(nxsupp/2+1:nxsupp,:,:) = w(nxsupp/2+1:nxsupp,:,:) + w(nxFp-nxsupp/2+1:nxFp,:,:);
    w(:,nysupp/2+1:nysupp,:) = w(:,nysupp/2+1:nysupp,:) + w(:,nyFp-nysupp/2+1:nyFp,:);
    
    % Return velocity field on original grid
    uout = u(nxsupp/2+1:nxFp-nxsupp/2,nysupp/2+1:nyFp-nysupp/2,:);
    vout = v(nxsupp/2+1:nxFp-nxsupp/2,nysupp/2+1:nyFp-nysupp/2,:);
    wout = w(nxsupp/2+1:nxFp-nxsupp/2,nysupp/2+1:nyFp-nysupp/2,:);