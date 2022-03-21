function E = getModelSpectrum(k,KE,L)
    % Inputs:
    %   nu --> molecular viscosity
    %   L --> Integral length scale
    %   C --> Model constant. Pope suggests C = 1.2, see pg. 233
    %   beta --> Model constant. Pope suggests beta = 5.2, see pg. 233
    %   KE --> kinetic energy
    
    % These use the high Re values given in Pope pg. 233, but these need to
    % be computed from the integral constraints of KE and dissipation
        cL = 6.78;
        ceta = 0.4;
    
    % use cL = 1 for KE formulation of spectrum
        cL = 1;
        
    % Nondimensional wavenumber
        kL = L*k;
        
    % Model coefficient chosen such that integrated spectrum equals total kinetic
    % energy
        C = 0.4843;
    
    % Energetic scales 
        fL = kL.^4./(cL + kL.^2).^(17/6); % Von Karman spectrum
        
    % Dissipative scales
%         eta = (nu^3/epsilon)^(1/4);
%         keta = kL.*eta/L;
%         feta = exp(-beta*((keta.^4 + ceta^4).^(1/4) - ceta));
        feta = ones(1,length(kL));
        disp('Warning: finite Re model spectrum is not implemented -- getModelSpectrum.m')
        
    % Model spectrum
        E = C*2*KE*L*fL.*feta;
