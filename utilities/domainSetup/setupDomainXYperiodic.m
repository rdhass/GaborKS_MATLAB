function [xLES,yLES,zLES,xQH,yQH,zQH,nxQH,nyQH,nzQH,dxQH,dyQH,dzQH,xF,yF,zF,xFp,yFp,nxsupp,nysupp,nzsupp,...
    kmin,kmax] = setupDomain(nxLES,nyLES,nzLES,Lx,Ly,Lz,nxLES_per_QH,nyLES_per_QH,nzLES_per_QH,nxF,nyF,nzF)
    % LES field
        dxLES = Lx/nxLES;
        dyLES = Ly/nyLES;
        dzLES = Lz/nzLES;

        xLES = 0:dxLES:Lx; 
        yLES = 0:dyLES:Ly; 
        zLES = dzLES/2:dzLES:Lz-dzLES/2;
    
    % QH mesh
        nxQH = nxLES/nxLES_per_QH;
        nyQH = nyLES/nyLES_per_QH;
        nzQH = nzLES/nzLES_per_QH;

        dxQH = Lx/nxQH;
        dyQH = Ly/nyQH;
        dzQH = Lz/nzQH;

        xQH = 0:dxQH:Lx;
        yQH = 0:dyQH:Ly;
        zQH = 0:dzQH:Lz; % Note that the QH mesh is defined from boundary to boundary
        
    % High resolution mesh
        dxF = Lx/nxF;
        dyF = Ly/nyF;
        dzF = Lz/nzF;
        
        xF = 0:dxF:Lx-dyF;
        yF = 0:dyF:Ly-dxF;
        zF = dzF/2:dzF:Lz-dzF/2;
        
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
        yen = xF(end)+nysupp/2*dyF;
        
        xFp = xst:dxF:xen;
        yFp = yst:dyF:yen;