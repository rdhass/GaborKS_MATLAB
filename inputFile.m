clc, clear all, close all

% Spatial domain
    Lx = 2*pi;
    Ly = 2*pi;
    Lz = 2*pi;

    % LES mesh
        nxLES = 16;
        nyLES = 16;
        nzLES = 16;

    % QH mesh
        nxLESperQH = 2;
        nyLESperQH = 2;
        nzLESperQH = 2;

    % High resolution mesh
        nxF = 4*nxLES;
        nyF = 4*nyLES;
        nzF = 4*nzLES;

% Large scale data 
    inputdir = '/Users/ryanhass/Documents/MATLAB/Lele_Research/Code development/GaborKS_V2_misc/data/';

% Where to write data
    outputdir = '/Users/ryanhass/Documents/MATLAB/Lele_Research/Code development/GaborKS_V2_misc/data/';

% Enrichment parameters
    nk = 10;
    ntheta = 10;
    kmin = 9;
    kmax = 32;
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