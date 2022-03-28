function [Uout,Vout,Wout] = extendVelocityToBoundaries(U,V,W,x,y,z,Lz)
    % Extend velocities to periodic boundary to make for easy interpolation
    % during mode evolution
    % Inputs:
    %   U,V,W --> large scale velocities rendered on LES mesh
    %   x,y,z --> 1D vectors specifying LES mesh plus boundary points
    
    nx = size(U,1);
    ny = size(U,2);
    nz = size(U,3);
    
    Uout = zeros(nx+1,ny+1,nz+2);
    Vout = zeros(nx+1,ny+1,nz+2);
    Wout = zeros(nx+1,ny+1,nz+2);
    
    Uout(1:nx,1:ny,2:nz+1) = U;
    Vout(1:nx,1:ny,2:nz+1) = V;
    Wout(1:nx,1:ny,2:nz+1) = W;
    
    % Enforce z-boundary conditions
    % Assume no-slip at bottom wall
    disp(['WARNING: Assuming no-slip velocity at solid wall for large scale velocity.',...
        ' This may lead to modeling errors.'])
    
    % Top wall is a slip wall
    %   w = 0; w is assumed to be an odd function
    %   dudz = 0 & dvdz = 0; u and v are assumed to be even functions
    % These imply the following
    %   u, v can be found by interpolation
    dz = z(3)-z(2);
    zvec = [0 dz/2:dz:Lz+dz/2];
    Uout(:,:,nz+2) = Uout(:,:,nz+1);
    Vout(:,:,nz+2) = Vout(:,:,nz+1);
    if size(x,1)~=size(U,1)
        x = transpose(x);
    end
    Uout = interp3(x,y,zvec,Uout,x,y,Lz,'spline');
    Vout = interp3(x,y,zvec,Vout,x,y,Lz,'spline');
    
    % Periodic BCs in x and y
    Uout(nx+1,:,:) = Uout(1,:,:);
    Uout(:,ny+1,:) = Uout(:,1,:);
    
    Vout(nx+1,:,:) = Vout(1,:,:);
    Vout(:,ny+1,:) = Vout(:,1,:);
    
    Wout(nx+1,:,:) = Wout(1,:,:);
    Wout(:,ny+1,:) = Wout(:,1,:);

   
    
