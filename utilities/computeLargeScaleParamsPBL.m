function [L,KE] = computeLargeScaleParamsPBL(U,V,W,Lx,Ly,Lz,x,y,z,cL)
    KE = 0.5*mean(U(:).^2 + V(:).^2 + W(:).^2);
    
    nx = size(U,1);
    ny = size(U,2);
    nz = size(U,3);
    
    vhat = fft(V,[],1)./nx;
    vhat = fft(vhat,[],2)./ny;
    R22 = ifft(ifft(vhat.*conj(vhat),[],1),[],2)*nx*ny;
    assert(max(abs(imag(R22)),[],'all') < 1e-14)
    
    xx = linspace(0,Lx/2,500)';
    yy = linspace(0,Ly/2,500);
   
    if size(x,1) ~= size(U,1)
        x = transpose(x);
    end
    
    Rup = interp3(x,y,z,R22,xx,yy,z,'spline');
    
    
    
