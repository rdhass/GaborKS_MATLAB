clc, clear all, close all
addpath('/Users/ryanhass/Documents/MATLAB/Lele_Research/Code development/GaborKS_V2/utilities')
load(['/Users/ryanhass/Documents/MATLAB/Lele_Research/Code development/',...
    'GaborKS_V2_misc/data/PBL/LargeScaleVelocity.mat'])

Lx = 2*pi;
Ly = pi;
Lz = 1;

nx = size(U,1);
ny = size(U,2);
nz = size(U,3);

dx = Lx/nx;
dy = Ly/ny;
dz = Lz/nz;

x = 0:dx:Lx-dx;
y = 0:dy:Ly-dy;
z = dz/2:dz:Lz-dz/2;

cL = 1;

[L,KE] = computeLargeScaleParamsPBL(U,V,W,Lx,Ly,x,y,z,cL);
figure
subplot(1,2,1)
plot(z,L,'-o')
subplot(1,2,2)
plot(z,KE,'-o')