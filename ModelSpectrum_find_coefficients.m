clc, clear all, close all

kL = linspace(0,1000,10000);

% Flow parameters
    nu = 0.0001;
    epsilon = 1;
    L = 1;
    KE = 1;
    eta = (nu^3/epsilon)^(1/4);
    
% Dissipation normalization for k
    keta = kL.*eta/L;

% Model parameters
    C = 1.2;
    beta = 5.2;
    
% These use the high Re values given in Pope pg. 233, but these need to
% be computed from the integral constraints of KE and dissipation
    ceta = 0.575;%0.4;
    cL = 1.758;%6.78;
    
% Energetic scales 
    fL = kL.^4./(cL + kL.^2).^(17/6); % Von Karman spectrum
    
% Dissipative scales
    feta = exp(-beta*((keta.^4 + ceta^4).^(1/4) - ceta));

% Model spectrum
    E = C*epsilon^(2/3)*L^(5/3)*fL.*feta;

% Compute epsilon and k from the spectrum
    kcheck = trapz(kL./L,E);
    epscheck = trapz(kL./L,2*nu*kL.^2./L^2.*E);

figure, hold on
plot(kL,E)
set(gca,'XScale','log','YScale','log')
%% Search for optimal cL and ceta
    cL = 10;
    ceta = 0;

iter = 0;
tol = 1e-1;
maxiter = 1000;
kcheck = 0;
epscheck = 0;
while kcheck < KE && epscheck < epsilon
    % Dissipative scales
        feta = exp(-beta*((keta.^4 + ceta^4).^(1/4) - ceta));
    while kcheck < KE && iter < maxiter
        cL = cL - 0.01;
        % Energetic scales 
            fL = kL.^4./(cL + kL.^2).^(17/6); % Von Karman spectrum

        % Model spectrum
            E = C*epsilon^(2/3)*L^(5/3)*fL.*feta;

        % Compute epsilon and k from the spectrum
            kcheck = trapz(kL./L,E);
            epscheck = trapz(kL./L,2*nu*kL.^2./L^2.*E);
            
        iter = iter + 1;
    end
    if iter == maxiter && kcheck < KE
        disp('Warning: computation for cL did not converge')
    elseif iter == maxiter
        disp('Warning: max iterations reached for cL, but kcheck < KE. There could be an error in the calculation')
    end
    iter = 0;
    while epscheck < epsilon && iter < maxiter
        ceta = ceta + 0.01;
        % Dissipative scales
            feta = exp(-beta*((keta.^4 + ceta^4).^(1/4) - ceta));
            
        % Model spectrum
            E = C*epsilon^(2/3)*L^(5/3)*fL.*feta;

        % Compute epsilon and k from the spectrum
            kcheck = trapz(kL./L,E);
            epscheck = trapz(kL./L,2*nu*kL.^2./L^2.*E);
        
        iter = iter + 1;
    end
    if iter == maxiter
        disp('Warning: computation for ceta did not converge')
    end
end
figure
subplot(2,1,1), hold on
plot(kL,fL)
plot(kL,feta)
set(gca,'YScale','log','XScale','log')
subplot(2,1,2)
plot(kL,E)
set(gca,'YScale','log','XScale','log')

