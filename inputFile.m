clc, clear all, close all

% Spatial domain
    Lx = 2*pi;
    Ly = 2*pi;
    Lz = 2*pi;

    % LES mesh
        nxLES = 32;
        nyLES = 32;
        nzLES = 32;

    % QH mesh
        nxLESperQH = 2;
        nyLESperQH = 2;
        nzLESperQH = 2;

    % High resolution mesh
        nxF = 8*nxLES;
        nyF = 8*nyLES;
        nzF = 8*nzLES;

% Large scale data 
    inputdir = '/work2/06632/ryanhass/stampede2/Enrichment/GaborKS_V2Data/HIT256_forced16/largeScales/';

% Where to write data
    outputdir = '/work2/06632/ryanhass/stampede2/Enrichment/GaborKS_V2Data/HIT256_forced16/';

% Enrichment parameters
    nk = 10;
    ntheta = 10;
    kmin = 11;
    kmax = 110;
    scalefact = sqrt(0.45); % A tuneable parameter to get the correct energy spectrum. 
                   % This is required due to the attenuation of mode energy 
                   % in the domain due to the window function
    ctau = 7;
    Anu = 2e-4; % SGS model constant
    epsilon = 1; % Large scale
    numolec = 0;
                   
% Misc parameters
    showGrid = false; % Graphically portray computational mesh? Good for debugging 
    doRenderIsotropic = true; % Render the pre-initialized isotropic modes?
    doRenderStrained = true; % Render the initialized anisotropic modes?
