function gradU = getLargeScaleGradient(U,V,W,Lx,Ly,Lz,xLES,yLES,zLES)
    % Compute large scale velocity gradient from large scale velocity
    % fields. Note, this requires knowledge of the boundary conditions
    % implemented in the LES solver
    % Inputs:
    %   U,V,W --> Large scale velocity
    %   Lx,Ly,Lz --> Domain size
    %   xLES,yLES,zLES --> 1D vectors defining LES grid (Note these vectors
    %                      include the boundary points!)
    %   z0 --> Roughness lengthscale used in simulation

    nx = size(U,1);
    ny = size(U,2);
    nz = size(U,3);
    
    dx = xLES(2)-xLES(1);
    dy = yLES(2)-yLES(1);
    dz = zLES(3)-zLES(2);
    
    gradU = zeros(3,3,nx,ny,nz);

    gradU(1,1,:,:,:) = ddx_hit(U,Lx,Ly,Lz);
    gradU(1,2,:,:,:) = ddy_hit(U,Lx,Ly,Lz);
    gradU(1,3,:,:,:) = ddzFD(U,dz);

    gradU(2,1,:,:,:) = ddx_hit(V,Lx,Ly,Lz);
    gradU(2,2,:,:,:) = ddy_hit(V,Lx,Ly,Lz);
    gradU(2,3,:,:,:) = ddzFD(V,dz);

    gradU(3,1,:,:,:) = ddx_hit(W,Lx,Ly,Lz);
    gradU(3,2,:,:,:) = ddy_hit(W,Lx,Ly,Lz);
    gradU(3,3,:,:,:) = ddzFD(W,dz);

    gradU2 = zeros(3,3,nx+1,ny+1,nz+2);
    gradU2(:,:,1:nx,1:ny,2:nz+1) = gradU;
    
    % dwdz can be estimated via 2nd order central difference and
    % noting that w is assumed to be an odd function at the top and bottom
    % walls
    gradU2(3,3,1:nx,1:ny,nz+2) = -2/dz*W(:,:,nz);
    gradU2(3,3,1:nx,1:ny,1) = 2/dz*W(:,:,1);
    
    % Assume u=v=0 at solid wall and use first order sided derivative.
    % Note, since the grid spacing at the wall is non-uniform (i.e. it goes
    % from dz/2 to dz), a second order stencil is invalid
    disp(['WARNING: assuming no-slip velocity at wall while computing large scale velocity gradient.',...
        ' This may lead to modeling errors.'])
    gradU2(1,3,1:nx,1:ny,1) = U(:,:,1)./dx;
    gradU2(2,3,1:nx,1:ny,1) = V(:,:,1)./dy;
    
    % Enforce periodic BCs in x and y
    gradU2(:,:,end,:,:) = gradU2(:,:,1,:,:);
    gradU2(:,:,:,end,:) = gradU2(:,:,:,1,:);
    
    % Enforce boundary conditions
    % Top wall is a slip wall
    %   w = 0; w is assumed to be an odd function
    %   dudz = 0 & dvdz = 0; u and v are assumed to be even functions
    % These imply the following
    %   dudx,dudy,dvdx,dvdy are assumed to be even functions and can be
    %                       found by interpolation 
    %   u, v can be found by interpolation
    %   dwdx, dwdy are assumed to be odd functions and can be found by
    %              interpolation
    %   dwdz can be computed from odd extension of w
    gradU2(1:2,1:2,:,:,nz+2) = gradU2(1:2,1:2,:,:,nz+1); % Even extension
    gradU2(3,1:2,:,:,nz+2) = -gradU2(3,1:2,:,:,nz+1); % Odd extension
    zvec = zeros(1,nz+1);
    zvec(1:nz) = zLES(2:nz+1);
    zvec(nz+1) = Lz;
    if size(gradU2,3) ~= size(xLES,1)
        xLES = transpose(xLES);
    end
    
    % Interpolate the values of dudx, dudy, dvdx, and dvdy using even extension
    for j = 1:2
        for i = 1:2
            gradU2(i,j,:,:,2:nz+2) = interp3(xLES,yLES,zvec,squeeze(gradU2(i,j,:,:,2:nz+2)),xLES,yLES,zLES(2:nz+2),'spline');
        end
    end
    
    % Interpolate the values of dwdx and dwdy using odd extension
    for j = 1:2
        gradU2(3,j,:,:,2:nz+2) = interp3(xLES,yLES,zvec,squeeze(gradU2(3,j,:,:,2:nz+2)),xLES,yLES,zLES(2:nz+2),'spline');
    end
    
    gradU = gradU2;
