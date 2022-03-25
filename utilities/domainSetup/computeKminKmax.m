function [kmin,kmax] = computeKminKmax(Lx,Ly,Lz,nxLES,nyLES,nzLES,nxF,nyF,nzF)

    kxNyqLES = getNyquist(Lx,nxLES);
    kyNyqLES = getNyquist(Ly,nyLES);
    kzNyqLES = getNyquist(Lz,nzLES);
    
    kxNyqF = getNyquist(Lx,nxF);
    kyNyqF = getNyquist(Ly,nyF);
    kzNyqF = getNyquist(Lz,nzF);
    
    kmin = min([kxNyqLES,kyNyqLES,kzNyqLES]);
    kmax = min([kxNyqF,kyNyqF,kzNyqF]);
