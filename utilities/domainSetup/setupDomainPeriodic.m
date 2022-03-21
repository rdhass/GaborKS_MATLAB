function [xLES,yLES,zLES,xQH,yQH,zQH,nxQH,nyQH,nzQH,dxQH,dyQH,dzQH,xF,yF,zF,xFp,yFp,zFp,nxsupp,nysupp,nzsupp] = setupDomainPeriodic(...
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

        xQH = 0:dxQH:Lx;
        yQH = 0:dyQH:Ly;
        zQH = 0:dzQH:Lz;
        
    % High resolution mesh
        dxF = Lx/nxF;
        dyF = Ly/nyF;
        dzF = Lz/nzF;
        
        xF = 0:dxF:Lx-dyF;
        yF = 0:dyF:Ly-dxF;
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
        
    % Do periodic extension of numerical mesh
        xst = -nxsupp/2*dxF;
        xen = xF(end)+nxsupp/2*dxF;
        
        yst = -nysupp/2*dyF;
        yen = xF(end)+nysupp/2*dyF;
        
        zst = -nzsupp/2*dzF;
        zen = zF(end)+nzsupp/2*dzF;
        
        xFp = xst:dxF:xen;
        yFp = yst:dyF:yen;
        zFp = zst:dzF:zen;