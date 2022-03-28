function gradU = getLargeScaleGradient(U,V,W,Lx,Ly,Lz)

    nx = size(U,1);
    ny = size(U,2);
    nz = size(U,3);
    
    gradU = zeros(3,3,nx,ny,nz);

    gradU(1,1,:,:,:) = ddx_hit(U,Lx,Ly,Lz);
    gradU(1,2,:,:,:) = ddy_hit(U,Lx,Ly,Lz);
    gradU(1,3,:,:,:) = ddz_hit(U,Lx,Ly,Lz);

    gradU(2,1,:,:,:) = ddx_hit(V,Lx,Ly,Lz);
    gradU(2,2,:,:,:) = ddy_hit(V,Lx,Ly,Lz);
    gradU(2,3,:,:,:) = ddz_hit(V,Lx,Ly,Lz);

    gradU(3,1,:,:,:) = ddx_hit(W,Lx,Ly,Lz);
    gradU(3,2,:,:,:) = ddy_hit(W,Lx,Ly,Lz);
    gradU(3,3,:,:,:) = ddz_hit(W,Lx,Ly,Lz);

    gradU2 = zeros(3,3,nx+1,ny+1,nz+1);
    gradU2(:,:,1:nx,1:ny,1:nz) = gradU;
    gradU2(:,:,end,:,:) = gradU2(:,:,1,:,:);
    gradU2(:,:,:,end,:) = gradU2(:,:,:,1,:);
    gradU2(:,:,:,:,end) = gradU2(:,:,:,:,1);

    gradU = gradU2;
