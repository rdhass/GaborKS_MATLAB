function [L,KE] = computeLargeScaleParams(U,V,W,Lx,Ly,xLES,yLES,zLES,nxLES,nyLES,zQHcent,nxQH,nyQH,nzQH,cL)
    % This computes the integral length scale and kinetic energy for the
    % turbulent half channel case and interpolates them to the center of
    % each QH region
    % Inputs:
    %   U, V, W --> Large scale velocity at each LES node
    %   Lx, Ly --> domain size in x and y
    %   xLES, yLES, zLES --> 1D vectors of LES node locations
    %   zQHcent --> 1D vector of QH center locations
    %   nxQH,nyQH,nzQH --> Number of QH regions
    %   cL --> tunable parameter for integral scale
    
    Uavg = mean(mean(U,1),2);
    Uavg = repmat(Uavg,[nxLES,nyLES,1]);
    U = U - Uavg;
    KE = 0.5*squeeze(mean(mean(U.^2 + V.^2 + W.^2,1),2));
    
    nx = size(U,1);
    ny = size(U,2);
    nz = size(U,3);
    
    vhat = fft(V,[],1)./nx;
    vhat = fft(vhat,[],2)./ny;
    R22 = ifft(ifft(vhat.*conj(vhat),[],1),[],2)*nx*ny;
    assert(max(abs(imag(R22)),[],'all') < 1e-14)
    
    % Interpolate R22 so we can find the 20% crossing 
    xx = linspace(0,Lx/2,500)';
    yy = linspace(0,Ly/2,500);
   
    if size(xLES,1) ~= size(U,1)
        xLES = transpose(xLES);
    end
    
    Rup = interp3(xLES(1:nx),yLES(1:ny),zLES(2:nz+1),R22,xx,yy,zLES(2:nz+1),'spline');
    
    
    Lxv = zeros(1,nz);
    Lyv = zeros(1,nz);
    for k = 1:nz
        f = squeeze(Rup(1,:,k)./Rup(1,1,k));
        g = squeeze(Rup(:,1,k)./Rup(1,1,k));
        yidx = find(f <= 0.2,1);
        xidx = find(g <= 0.2,1);
        
        Lyv(k) = trapz(yy(1:yidx),f(1:yidx));
        Lxv(k) = trapz(xx(1:xidx),g(1:xidx));
    end
    
    L = cL*sqrt(Lyv.*Lxv);
    
    % Interpolate L and KE to QH centers
    L = interp1(zLES(2:nz+1),L,zQHcent,'spline');
    KE = interp1(zLES(2:nz+1),KE,zQHcent,'spline');
    
    % Replicate vectors to appropriate sized arrays
    L = reshape(L,[1,1,nzQH]);
    KE = reshape(KE,[1,1,nzQH]);
    
    L = repmat(L,[nxQH,nyQH,1]);
    KE = repmat(KE,[nxQH,nyQH,1]);
    
    
