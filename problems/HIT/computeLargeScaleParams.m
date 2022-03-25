function [L11,KE] = computeLargeScaleParamsHIT(U,V,W,Lx,Ly,Lz)
    [E,kline] = get_3Denergy_spectrum(U,V,W,Lx,Ly,Lz,size(U,1));
    KE = 0.5*mean(U(:).^2 + V(:).^2 + W(:).^2);
    L11 = 3*pi/(4*KE)*trapz(kline(2:end),E(2:end)./kline(2:end)'); % Pope eqn (6.260)