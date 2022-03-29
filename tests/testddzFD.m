clc, clear all, close all

nz = 100;
nx = 100;
ny = 100;

Lx = 1;
Ly = 1;
Lz = 1;

dx = Lx/nx;
dy = Ly/ny;
dz = Lz/nz;

x = 0:dx:Lx-dx;
y = 0:dy:Ly-dy;
z = dz/2:dz:Lz-dz/2;

[x,y,z] = ndgrid(x,y,z);

f = sin(2*pi*z);
fpTrue = 2*pi*cos(2*pi*z);
fpNum = ddzFD(f,dz);

figure, hold on
plot(squeeze(z(1,1,:)),squeeze(f(1,1,:)),'DisplayName','$f(z)$')
plot(squeeze(z(1,1,:)),squeeze(fpNum(1,1,:)),'DisplayName',"$f'_{num}$")
plot(squeeze(z(1,1,:)),squeeze(fpTrue(1,1,:)),'--','DisplayName',"$f'_{exact}$")
legend('Interpreter','Latex')