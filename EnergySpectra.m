clc, clear all, close all
addpath('./functions')
addpath('/Users/ryanhass/Documents/MATLAB/Lele_Research/matlabFunctions')
load('LargeScaleVelocity.mat')

H = 2*pi;
[E,kline] = get_3Denergy_spectrum(U,V,W,H,H,H,16);
figure, hold on

xlim([8e-1 5e1])
ylim([1e-5 1e0])
KE = trapz(kline,E);
L11 = 3*pi/(4*KE)*trapz(kline(2:end),E(2:end)./kline(2:end)');

% Model spectrum
kmod = linspace(0.1,32,1000);
Cvk = 0.4843;
Emod = Cvk.*kmod.^4./(1+kmod.^2).^(17/6);
Emod2 = getModelSpectrum(kmod,KE,L11);

% Gabor-induced field (isotropic)
nthetavec = [10 20];
nkvec = [10 20 30];
Efact = 1;
pid = 0;
for ntid = 1:length(nthetavec)
    ntheta = nthetavec(ntid);
    for nkid = 1:length(nkvec)
        pid = pid + 1;
        subplot(2,3,pid), hold on
        plot(kline,E,'DisplayName','Large scale field')
        plot(kmod,Emod,'DisplayName','Model spectrum')
        plot(kmod,Emod2/(2*KE*L11),'DisplayName','Model spectrum 2')

        nk = nkvec(nkid);
        load(['IsotropicVelocity_ntheta',num2str(ntheta),'_nk',num2str(nk),'.mat'])
        [EG,kG] = get_3Denergy_spectrum(u,v,w,H,H,H,64);
        plot(kG,Efact*EG./(2*KE*L11),'DisplayName','Gabor-induced field (isotropic modes)')
        
        if nk == 10 && ntheta == 10
            
        end

        % Compute spectrum of Gabor modes based on mode values instead of induced
        % velocity
        load(['GaborModes_ntheta',num2str(ntheta),'_nk',num2str(nk),'.mat'])
        nmodes = length(gmxloc);
        kmin = 9;
        kmax = 32;

        nQH = 8^3;
        kedge = logspace(log10(kmin),log10(kmax),nk+1);
        kmag = (kedge(1:end-1) + kedge(2:end))/2;
        EG2 = zeros(1,length(kedge)-1);
        for n = 1:nmodes
            k = sqrt(kx(n)^2 + ky(n)^2 + kz(n)^2);
            bin = find(kedge >= k,1)-1;
            EG2(bin) = EG2(bin) + 0.5*(uhatR(n)^2 + uhatI(n)^2 + ...
                                       vhatR(n)^2 + vhatI(n)^2 + ...
                                       whatR(n)^2 + whatI(n)^2);
        end
        for kid = 1:length(kmag)
            EG2(kid) = EG2(kid)/(kedge(kid+1) - kedge(kid));
        end
        EG2 = EG2./nQH;
        plot(kmag,EG2./(2*KE*L11),'.')
        set(gca,'YScale','log','XScale','log')
        xlim([8e-1 1e2])
        ylim([1e-3 1e0])
    end
end
subplot(2,3,1)
legend('Interpreter','Latex')

ntheta = 10;
nk = 10;
load(['StrainedVelocity_ntheta',num2str(ntheta),'_nk',num2str(nk),'.mat'])
[EGstr,~] = get_3Denergy_spectrum(u,v,w,H,H,H,64);
plot(kG,Efact*EGstr./(2*KE*L11),'DisplayName','Gabor-induced field (strained modes)')