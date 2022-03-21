function [u,v,w] = renderVelocity(uhatR,uhatI,vhatR,vhatI,whatR,whatI,kx,ky,kz,gmx,gmy,gmz,nxsupp,nysupp,nzsupp,xF,yF,zF)
    % Inputs:
    %   uhatR, uhatI, ... --> complex velocity amplitudes associated with
    %                         each Gabor mode. Array size: nmodesX1
    %   kx, ky, kz --> wave-vector components for each GM. Array size: nmodesX1
    %   gmx, gmy, gmz --> Physical location of each GM. Array size: nmodesX1
    %   nxsupp, nysupp, nzsupp --> Number of fine grid points (on each side of a given GM) 
    %                              spanned by the window function. Scalar
    %                              (integer)
    %   xF, yF, zF --> Eulerian grid for fine mesh. Array size: 1XnF
    % Outputs:
    %   u, v, w --> Physical space velocities. Array size: nxF X nyF X nzF
        
    dxF = xF(2)-xF(1);
    dyF = yF(2)-yF(1);
    dzF = zF(2)-zF(1);
    
    nxF = length(xF);
    nyF = length(yF);
    nzF = length(zF);
    
    u = zeros(nxF,nyF,nzF);
    v = u;
    w = u;
    
    [xF,yF,zF] = ndgrid(xF,yF,zF);
    
    wxSupport = (nxsupp + 1)*dxF;
    wySupport = (nysupp + 1)*dyF;
    wzSupport = (nzsupp + 1)*dzF;
    
    nmodes = length(kx);
 
    % Given a mode location and its support width we can convert this
    % information to indices on the fine grid 
    for n = 1:nmodes
        ist = max([1,ceil(gmx(n)/dxF)  - nxsupp/2]);
        ien = min([nxF,floor(gmx(n)/dxF) + nxsupp/2]);
        
        jst = max([1,ceil(gmy(n)/dyF)  - nysupp/2]);
        jen = min([nyF,floor(gmy(n)/dyF) + nysupp/2]);
        
        kst = max([1,ceil(gmz(n)/dzF)  - nzsupp/2]);
        ken = min([nzF,floor(gmz(n)/dzF) + nzsupp/2]);
        
        x = xF(ist:ien,jst:jen,kst:ken);
        y = yF(ist:ien,jst:jen,kst:ken);
        z = zF(ist:ien,jst:jen,kst:ken);
        
        kdotx = kx(n).*(x-gmx(n)) + ky(n).*(y-gmy(n)) + kz(n).*(z-gmz(n));
        cs = cos(kdotx);
        ss = sin(kdotx);
        
        % Window function
        fx = cos(pi*(x - gmx(n))/wxSupport);
        fy = cos(pi*(y - gmy(n))/wySupport);
        fz = cos(pi*(z - gmz(n))/wzSupport);
        f = fx.*fy.*fz;
        
        % Resulting velocity field
        u(ist:ien,jst:jen,kst:ken) = u(ist:ien,jst:jen,kst:ken) + f.*(2*uhatR(n).*cs - 2*uhatI(n).*ss);
        v(ist:ien,jst:jen,kst:ken) = v(ist:ien,jst:jen,kst:ken) + f.*(2*vhatR(n).*cs - 2*vhatI(n).*ss);
        w(ist:ien,jst:jen,kst:ken) = w(ist:ien,jst:jen,kst:ken) + f.*(2*whatR(n).*cs - 2*whatI(n).*ss);
    end