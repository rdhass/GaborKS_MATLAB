function [L,KE] = computeLargeScaleParamsPBL(U,V,W,Lx,Ly,x,y,z,cL)
    % This computes the integral length scale and kinetic energy for the
    % turbulent half channel case
    % Inputs:
    %   U, V, W --> Large scale velocity at each LES node
    %   Lx, Ly --> domain size in x and y
    %   x, y, z --> 1D vectors of LES node locations
    %   cL --> tunable parameter for integral scale
    
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
   
    if size(x,1) ~= size(U,1)
        x = transpose(x);
    end
    
    Rup = interp3(x,y,z,R22,xx,yy,z,'spline');
    
    
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
    
    
