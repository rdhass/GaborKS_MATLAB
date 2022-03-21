clc, clear all, close all

% Load data from a 256^3 LES whos dissipation is set exclusively by the SGS
% model, i.e. inifinite Re simulation. 
load('E_L11_dissipation_and_KE_t015627.mat')

% Is it better to specify the large scales spectrum via its energy
% dissipation rate or its kinetic energy?

% Plot computed spectrum which is normalized to give KE
figure, hold on
plot(kline,E,'LineWidth',1.3,'DisplayName','Computed from $256^3$ LES');

% Plot two model spectra:
%   1. Normalized to give correct KE
%   2. The one reported in Aditya's thesis

%   1:
    Cvk = 0.4843;
    kL = kline*L11;
    E1 = Cvk*KE*2*L11*kL.^4./(1 + kL.^2).^(17/6);
    kLL = kline*L11L;
    E1L = Cvk*KEL*2*L11*kLL.^4./(1 + kLL.^2).^(17/6);
    plot(kline,E1,'r','LineWidth',1.3,'DisplayName','$E_1(k)$ from full field')
    plot(kline,E1L,'r--','LineWidth',1.3,'DisplayName','$E_1(k)$ from scale-speprated large scales')
    
%   2:
    alpha = 1.7;
    E2 = alpha*epsilon^(2/3)*L11^(-5/3)*kL.^4./(1 + kL.^2).^(17/6);
    E2L = alpha*epsilon^(2/3)*L11L^(-5/3)*kLL.^4./(1 + kLL.^2).^(17/6);
    plot(kline,E2,'k','LineWidth',1.3,'DisplayName','$E_2(k)$ from full field')
    plot(kline,E2L,'k--','LineWidth',1.3,'DisplayName','$E_2(k)$ from scale-speprated large scales')
    
kcoLES = (2/3)*16;
kcoLS = (2/3)*8;
plot([kcoLES,kcoLES],[1e-6 1e-1],'k','DisplayName','$k_{co}$ corresponding to $32^3$ LES')
plot([kcoLS,kcoLS],[1e-6 1e-1],'k--','DisplayName','$k_{co}$ corresponding to scale-separated field')
    
set(gca,'YScale','log','XScale','log')
ylim([1e-6 1e-1])
xlim([8e-1 3e2])
legend('Interpreter','Latex')
xlabel('$k$','Interpreter','Latex')
ylabel('$E(k)$','Interpreter','Latex')