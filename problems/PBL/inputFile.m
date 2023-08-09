clc, clear all, close all

% Spatial domain
    Lx = 2*pi;
    Ly = pi;
    Lz = 1;

    % LES mesh
        nxLES = 192;
        nyLES = 192;
        nzLES = 64;

    % QH mesh
        nxLESperQH = 2;
        nyLESperQH = 2;
        nzLESperQH = 2;

    % High resolution mesh
        nxF = 2*nxLES;
        nyF = 2*nyLES;
        nzF = 2*nzLES;

% Large scale data 
    inputdir = '/work2/06632/ryanhass/stampede2/Enrichment/GaborKS_V2Data/PBL64/';

% Where to write data
    outputdir = '/work2/06632/ryanhass/stampede2/Enrichment/GaborKS_V2Data/PBL64/';

% Enrichment parameters
    nk = 10;
    ntheta = 10;
    scalefact = sqrt(0.45); % A tuneable parameter to get the correct energy spectrum. 
                   % This is required due to the attenuation of mode energy 
                   % in the domain due to the window function
    ctau = 7;
    Anu = 1e-4; % SGS model constant
    epsilon = 1; % Large scale
    numolec = 0;
    cL = 1; % Parameter to adjust the integral length scale of the large scale field
                   
% Misc parameters
    showGrid = false; % Graphically portray computational mesh? Good for debugging 
    doRenderIsotropic = true; % Render the pre-initialized isotropic modes?
    doRenderStrained = true; % Render the initialized anisotropic modes?
