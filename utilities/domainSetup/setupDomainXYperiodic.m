function [xLES,yLES,zLES,xQHedge,yQHedge,zQHedge,xQHcent,yQHcent,zQHcent,nxQH,nyQH,nzQH,dxQH,dyQH,dzQH,...
    xF,yF,zF,xFp,yFp,nxsupp,nysupp,nzsupp,kmin,kmax] = setupDomainXYperiodic(...
    nxLES,nyLES,nzLES,Lx,Ly,Lz,nxLES_per_QH,nyLES_per_QH,nzLES_per_QH,nxF,nyF,nzF)
    % This routine sets up the spatial mesh. Note that the "LES" and "QH"
    % grids include the boundaries which is different than the LES solver
    % which excludes one side of periodic boundaries and excludes the top
    % and bottom wall for wall-bounded flows
    
    % LES field
        dxLES = Lx/nxLES;
        dyLES = Ly/nyLES;
        dzLES = Lz/nzLES;

        xLES = 0:dxLES:Lx; 
        yLES = 0:dyLES:Ly; 
        zLES = [0 dzLES/2:dzLES:Lz-dzLES/2 Lz];
    
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