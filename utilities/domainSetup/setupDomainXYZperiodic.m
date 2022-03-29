function [xLES,yLES,zLES,xQHedge,yQHedge,zQHedge,xQHcent,yQHcent,zQHcent,nxQH,nyQH,nzQH,dxQH,dyQH,dzQH,...
    xF,yF,zF,xFp,yFp,zFp,nxsupp,nysupp,nzsupp, kmin,kmax] = setupDomainPeriodic(...
    nxLES,nyLES,nzLES,Lx,Ly,Lz,nxLES_per_QH,nyLES_per_QH,nzLES_per_QH,nxF,nyF,nzF)
    % LES field
        dxLES = Lx/nxLES;
        dyLES = Ly/nyLES;
        dzLES = Lz/nzLES;

        xLES = 0:dxLES:Lx; 
        yLES = 0:dyLES:Ly; 
        zLES = 0:dzLES:Lz;
    
    % QH mesh
        nxQH = nxLES/nxLES_per_QH;
        nyQH = nyLES/nyLES_per_QH;
        nzQH = nzLES/nzLES_per_QH;

        dxQH = Lx/nxQH;
        dyQH = Ly/nyQH;
        dzQH = Lz/nzQH;

        xQHedge = 0:dxQH:Lx;
        yQHedge = 0:dyQH:Ly;
        zQHedge = 0:dzQH:Lz;
        
        xQHcent = 0.5*(xQHedge(1:end-1) + xQHedge(2:end));
        yQHcent = 0.5*(yQHedge(1:end-1) + yQHedge(2:end));
        zQHcent = 0.5*(zQHedge(1:end-1) + zQHedge(2:end));
        
    % High resolution mesh
        dxF = Lx/nxF;
        dyF = Ly/nyF;
        dzF = Lz/nzF;
        
        xF = 0:dxF:Lx-dxF;
        yF = 0:dyF:Ly-dyF;
        zF = 0:dzF:Lz-dzF;
        
    % Domain of influence of each mode
        nxF_per_QH = nxF/nxQH;
        nyF_per_QH = nyF/nyQH;
        nzF_per_QH = nzF/nzQH;
        
        % The problem is formulated such that the domain of influence of
        % each mode is 2 QH regions in extent.
        nxsupp = 2*nxF_per_QH;
        nysupp = 2*nyF_per_QH;
        nzsupp = 2*nzF_per_QH;
        
    % Compute kmin and kmax
        [kmin,kmax] = computeKminKmax(Lx,Ly,Lz,nxLES,nyLES,nzLES,nxF,nyF,nzF);
        
    % Do periodic extension of numerical mesh
        xst = -nxsupp/2*dxF;
        xen = xF(end)+nxsupp/2*dxF;
        
        yst = -nysupp/2*dyF;
        yen = yF(end)+nysupp/2*dyF;
        
        zst = -nzsupp/2*dzF;
        zen = zF(end)+nzsupp/2*dzF;
        
        xFp = xst:dxF:xen;
        yFp = yst:dyF:yen;
        zFp = zst:dzF:zen;