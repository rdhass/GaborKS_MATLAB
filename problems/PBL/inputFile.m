clc, clear all, close all

% Spatial domain
    Lx = 2*pi;
    Ly = pi;
    Lz = 1;

    % LES mesh
        nxLES = 32;
        nyLES = 32;
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
    inputdir = '/Users/ryanhass/Documents/MATLAB/Lele_Research/Code development/GaborKS_V2_misc/data/PBL/';

% Where to write data
    outputdir = '/Users/ryanhass/Documents/MATLAB/Lele_Research/Code development/GaborKS_V2_misc/data/PBL/';

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
    showGrid = true; % Graphically portray computational mesh? Good for debugging 
    doRenderIsotropic = true; % Render the pre-initialized isotropic modes?
    doRenderStrained = true; % Render the initialized anisotropic modes?