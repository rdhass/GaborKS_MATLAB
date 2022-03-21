clc, clear all, close all
addpath('./data')
load('IsotropicVelocity_ntheta10_nk10.mat')
H = 2*pi;
h = H/nxF;
x = (0:h:H-h)';
y = 0:h:H-h;
x = x*ones(1,nyF);
y = ones(nxF,1)*y;

figure
subplot(2,3,1)
pcolor(x,y,squeeze(u(:,:,end/2)))
colorbar
shading interp
daspect([1 1 1])

subplot(2,3,2)
pcolor(x,y,squeeze(v(:,:,end/2)))
colorbar
shading interp
daspect([1 1 1])

subplot(2,3,3)
pcolor(x,y,squeeze(w(:,:,end/2)))
colorbar
shading interp
daspect([1 1 1])

load('StrainedVelocity_ntheta10_nk10.mat')
subplot(2,3,4)
pcolor(x,y,squeeze(u(:,:,end/2)))
colorbar
shading flat
daspect([1 1 1])

subplot(2,3,5)
pcolor(x,y,squeeze(v(:,:,end/2)))
colorbar
shading flat
daspect([1 1 1])

subplot(2,3,6)
pcolor(x,y,squeeze(w(:,:,end/2)))
colorbar
shading flat
daspect([1 1 1])
%%
% Plot mode locations
load('GaborModes.mat')

figure, hold on
plot(gmxloc,gmyloc,'.')