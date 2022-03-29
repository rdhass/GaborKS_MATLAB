clc, clear all, close all
addpath('../utilities')
addpath('../utilities/derivatives')
addpath('/Users/ryanhass/Documents/MATLAB/Lele_Research/matlabFunctions')
nz = 200;
nx = 100;
ny = 100;

Lx = 2*pi;
Ly = pi;
Lz = 1;

dx = Lx/nx;
dy = Ly/ny;
dz = Lz/nz;

x = 0:dx:Lx-dx;
y = 0:dy:Ly-dy;
z = 0:dz:Lz;

[x,y,z] = ndgrid(x,y,z);

u = (sin(2*x) + 0.2*cos(12*x)).*sin(2*y).*(z.^3 + 2*z.^2 - 0.5*z - 3);
v = (cos(3*x) + 0.1*sin(15*x)).*cos(3*y).*(0.2*z.^3 + z.^2 - 2*z);
w = cos(3*x).*cos(2*y).*cos(z);

figure
subplot(1,2,1)
pcolor(squeeze(w(:,end/2,:))')
shading flat
colorbar

[~,~,wout] = enforceNoPenBCxyPeriodic(u,v,w,Lx,Ly,Lz,squeeze(z(1,1,:)));
subplot(1,2,2)
pcolor(squeeze(wout(:,end/2,:))')
shading flat
colorbar

figure
for j = 1:ny
    plot(squeeze(x(:,j,1)),squeeze(w(:,j,1))), hold on
    plot(squeeze(x(:,j,1)),squeeze(wout(:,j,1))), hold off
    grid on
    ylim([-0.1 0.1])
    frame=getframe;
end