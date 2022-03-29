clc, clear all, close all
addpath('../')

% Set search paths
    searchPaths

% Read input file to set problem-specific parameters. 
    inputFile

% Setup the spatial domain
    [xLES,yLES,zLES,...                                                                 % LES grid
     xQHedge,yQHedge,zQHedge,xQHcent,yQHcent,zQHcent,nxQH,nyQH,nzQH,dxQH,dyQH,dzQH,...  % QH grid
     xF,yF,zF,xFp,yFp,zFp,...                                                           % High res grid and its periodic extension
     nxsupp,nysupp,nzsupp,kmin,kmax] ...                                                % Mode support
     = setupDomainXYZperiodic(nxLES,nyLES,nzLES,Lx,Ly,Lz,nxLESperQH,...
            nyLESperQH,nzLESperQH,nxF,nyF,nzF);
    
% Read in large scale flow field and compute integral scales
    load([inputdir,'LargeScaleVelocity.mat'])
    gradU = getLargeScaleGradient(U,V,W,Lx,Ly,Lz);
    [L,KE] = computeLargeScaleParams(U,V,W,Lx,Ly,Lz,nxQH,nyQH,nzQH);
    [U,V,W] = extendVelocityToBoundaries(U,V,W);
    
% Generate isotropic modes
    % First assign global variables
        global uhatR vhatR whatR;
        global uhatI vhatI whatI;
        global kx ky kz;
        global gmxloc gmyloc gmzloc;
        
    [uhatR,uhatI,vhatR,vhatI,whatR,whatI,kx,ky,kz,gmxloc,gmyloc,gmzloc] ...
     = initializeIsotropicModes(nxQH,nyQH,nzQH,xQHedge,yQHedge,zQHedge,dxQH,dyQH,dzQH,...
       nk,ntheta,kmin,kmax,scalefact,KE,L);
    disp('Finished initializing isotropic Gabor modes. Writing data to disk')
    save([outputdir,'GaborModes_ntheta',num2str(ntheta),'_nk',num2str(nk),'_V2.mat'],...
        'uhatR','uhatI','vhatR','vhatI','whatR','whatI','kx','ky','kz',...
        'gmxloc','gmyloc','gmzloc','nxQH','nyQH','nzQH','nxLES','nyLES','nzLES','nxF','nyF','nzF',...
        'nxsupp','nysupp','nzsupp','xF','yF','zF','xFp','yFp','zFp')
    
% Visualize computational mesh with Gabor mode locations shown in one QH
% region
    if showGrid
        plotGrid(xLES,yLES,zLES,nxQH,nyQH,nzQH,xQHedge,yQHedge,zQHedge,xF,yF,zF,...
            xFp,yFp,gmxloc,gmyloc,gmzloc,nxsupp,nysupp)
    end
%%
% Render velocity field
    if doRenderIsotropic
        disp('Rendering isotropic velocity field')
        [u,v,w] = renderVelocityXYZperiodic(uhatR,uhatI,vhatR,vhatI,whatR,whatI,...
            kx,ky,kz,gmxloc,gmyloc,gmzloc,nxsupp,nysupp,nzsupp,xFp,yFp,zFp);
        disp('Velocity rendered. Saving data to disk')
        save([outputdir,'IsotropicVelocity_ntheta',num2str(ntheta),'_nk',num2str(nk),'_V2.mat'],...
            'u','v','w','nxF','nyF','nzF')
    end
    %%
% Strain the modes
    strainModes(xLES,yLES,zLES,xQHcent,yQHcent,zQHcent,gradU,U,V,W,...
        gmxloc,gmyloc,gmzloc,kx,ky,kz,uhatR,uhatI,vhatR,vhatI,whatR,whatI,...
        Anu,KE,L,numolec,ctau);
    disp('Finished evolving Gabor modes. Writing data to disk')
    save([outputdir,'GaborModesAfterStraining_ntheta',num2str(ntheta),'_nk',num2str(nk),'_V2.mat'],...
        'uhatR','uhatI','vhatR','vhatI','whatR','whatI','kx','ky','kz',...
    'gmxloc','gmyloc','gmzloc','nxQH','nyQH','nzQH','nxLES','nyLES','nzLES','nxF','nyF','nzF',...
    'nxsupp','nysupp','nzsupp','xF','yF','zF')
%%
% Render velocity field
    if doRenderStrained
        disp('Begin rendering velocity')
        [u,v,w] = renderVelocityXYZperiodic(uhatR,uhatI,vhatR,vhatI,whatR,whatI,kx,ky,kz,...
            gmxloc,gmyloc,gmzloc,nxsupp,nysupp,nzsupp,xFp,yFp,zFp);
        disp('Writing velocity to disk')
        save([outputdir,'StrainedVelocity_ntheta',num2str(ntheta),'_nk',num2str(nk),'_V2.mat'],...
            'u','v','w','nxF','nyF','nzF')
    end
