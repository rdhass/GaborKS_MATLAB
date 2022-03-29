function [uout,vout,wout] = enforceNoPenBCxyPeriodic(u,v,w,Lx,Ly,Lz,z)
    % Superpose irrotational velocity field on Gabor-induced field that
    % exactly cancels the non-zero penetration velocity at the wall
    % Inputs:
    %   u,v,w --> Gabor-induced velocity that does not satisfy no-pen BC
    %   Lx,Ly,Lz --> domain size
    %   z --> 1D vector defining z mesh points. IMPORTANT: make sure this
    %         is defined to include domain boundaries, i.e. z = 0:dz:Lz.
    %         and make sure this coincides with how u, v, and w are defined
    
    nx = size(u,1);
    ny = size(u,2);
    nz = size(u,3)-1;
    
    assert(length(z) == nz+1);
    dz = z(2)-z(1);
    
    kx = fftshift(-nx/2:nx/2-1)'*2*pi./Lx;
    ky = fftshift(-ny/2:ny/2-1)*2*pi./Ly;
    kx = kx*ones(1,ny);
    ky = ones(nx,1)*ky;
    
    k2D = sqrt(kx.^2 + ky.^2);
    
    what = fft(fft(w,[],1)./nx,[],2)./ny;
    
    phiHat = zeros(nx,ny,nz+1);
    what0 = squeeze(what(:,:,1));
    what0 = repmat(what0,[1,1,nz+1]);
    k2D = repmat(k2D,[1,1,nz+1]);
    z = reshape(z,1,1,nz+1);
    z = repmat(z,[nx,ny,1]);
    phiHat = -what0./k2D.*exp(-k2D.*z);
    phiHat(1,1,:) = 0;
    
    phi = ifft(ifft(phiHat,[],2,'symmetric')*ny,[],1,'symmetric')*nx;
%     assert(false,'Need to figure out why irrotational field does not satisfy boundary conditions')
    uout = u - ddx_hit(phi,Lx,Ly,Lz);
    vout = v - ddy_hit(phi,Lx,Ly,Lz);
    wout = w - ddzFD(phi,dz);