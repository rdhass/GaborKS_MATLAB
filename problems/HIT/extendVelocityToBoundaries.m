function [Uout,Vout,Wout] = extendVelocityToBoundaries(U,V,W)
    % Extend velocities to periodic boundary to make for easy interpolation
    % during mode evolution
    
    nx = size(U,1);
    ny = size(U,2);
    nz = size(U,3);
    
    Uout = zeros(nx+1,ny+1,nz+1);
    Vout = zeros(nx+1,ny+1,nz+1);
    Wout = zeros(nx+1,ny+1,nz+1);
    
    Uout(1:nx,1:ny,1:nz) = U;
    Uout(nx+1,:,:) = Uout(1,:,:);
    Uout(:,ny+1,:) = Uout(:,1,:);
    Uout(:,:,nz+1) = Uout(:,:,1);
    
    Vout(1:nx,1:ny,1:nz) = V;
    Vout(nx+1,:,:) = Vout(1,:,:);
    Vout(:,ny+1,:) = Vout(:,1,:);
    Vout(:,:,nz+1) = Vout(:,:,1);
    
    Wout(1:nx,1:ny,1:nz) = W;
    Wout(nx+1,:,:) = Wout(1,:,:);
    Wout(:,ny+1,:) = Wout(:,1,:);
    Wout(:,:,nz+1) = Wout(:,:,1);
