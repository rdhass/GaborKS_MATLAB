clc, clear all, close all
addpath('./data')
addpath('./functions')

% The large scale velocity field is scale-separated, meaning low-pass
% filtered LES data rendered on the original LES mesh (as opposed to a
% downsampled version capable of resolving the filtered field).
load('LargeScaleVelocity.mat')

% Domain setup
    H = 2*pi;
    nxLES = 16;
    nyLES = 16;
    nzLES = 16;

gradU = zeros(3,3,nxLES,nyLES,nzLES);

gradU(1,1,:,:,:) = ddx_hit(U,H,H,H);
gradU(1,2,:,:,:) = ddy_hit(U,H,H,H);
gradU(1,3,:,:,:) = ddz_hit(U,H,H,H);

gradU(2,1,:,:,:) = ddx_hit(V,H,H,H);
gradU(2,2,:,:,:) = ddy_hit(V,H,H,H);
gradU(2,3,:,:,:) = ddz_hit(V,H,H,H);

gradU(3,1,:,:,:) = ddx_hit(W,H,H,H);
gradU(3,2,:,:,:) = ddy_hit(W,H,H,H);
gradU(3,3,:,:,:) = ddz_hit(W,H,H,H);

gradU2 = zeros(3,3,nxLES+1,nyLES+1,nzLES+1);
gradU2(:,:,1:nxLES,1:nyLES,1:nzLES) = gradU;
gradU2(:,:,end,:,:) = gradU2(:,:,1,:,:);
gradU2(:,:,:,end,:) = gradU2(:,:,:,1,:);
gradU2(:,:,:,:,end) = gradU2(:,:,:,:,1);

gradU = gradU2;

save('./data/LargeScaleGradient.mat','gradU')
