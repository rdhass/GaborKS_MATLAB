function plotGrid(xLES,yLES,nxQH,nyQH,nzQH,xQH,yQH,xF,yF,xFp,yFp,gmxloc,gmyloc,nxsupp,nysupp)
    nxLES = length(xLES);
    nyLES = length(yLES);
    xLES2D = xLES'*ones(1,nyLES);
    yLES2D = ones(nxLES,1)*yLES;

    nxQH = length(xQH)-1;
    nyQH = length(yQH)-1;
    xQH2D = xQH'*ones(1,nyQH+1);
    yQH2D = ones(nxQH+1,1)*yQH;
    
    nxF = length(xF);
    nyF = length(yF);
    xF2D = xF'*ones(1,nyF);
    yF2D = ones(nxF,1)*yF;
    dxF = xF(2)-xF(1);
    dyF = yF(2)-yF(1);
    
    nxFp = length(xFp);
    nyFp = length(yFp);
    xFp2D = xFp'*ones(1,nyFp);
    yFp2D = ones(nxFp,1)*yFp;
    
    gmxloc = reshape(gmxloc,nxQH,nyQH,nzQH,[]);
    gmyloc = reshape(gmyloc,nxQH,nyQH,nzQH,[]);
     
    figure
    subplot(1,2,1), hold on
    plot(xLES2D(1,1),yLES2D(1,1),'rx','DisplayName','LES nodes')
    plot(xQH2D(1,1),yQH2D(1,1),'bo','DisplayName','QH boundary')
    plot(xF2D(1,1),yF2D(1,1),'g.','DisplayName','Fine grid')
    plot(squeeze(gmxloc(4,4,4,1)),squeeze(gmyloc(4,4,4,1)),'k.','DisplayName','Gabor mode locations')
    plot(xLES2D,yLES2D,'rx','HandleVisibility','off')
    plot(xQH2D,yQH2D,'bo','HandleVisibility','off')
    plot(xF2D,yF2D,'g.','HandleVisibility','off')
    plot(squeeze(gmxloc(4,4,4,:)),squeeze(gmyloc(4,4,4,:)),'k.','HandleVisibility','off')
    xticks([0 pi/2 pi 3*pi/2 2*pi])
    yticks([0 pi/2 pi 3*pi/2 2*pi])
    xticklabels({'0', '\pi/2','\pi','3\pi/2','2\pi'})
    yticklabels({'0', '\pi/2','\pi','3\pi/2','2\pi'})
    grid on
    legend('Interpreter','Latex')
    daspect([1 1 1])
    
    subplot(1,2,2), hold on
    mrksz1 = 5;
    mrksz2 = 3;
    plot(xF2D(1,1),yF2D(1,1),'k.','MarkerSize',mrksz1,'DisplayName','Fine grid')
    plot(xF2D,yF2D,'k.','MarkerSize',mrksz1,'HandleVisibility','off')
    plot(xFp2D(1,1),yFp2D(1,1),'b.','MarkerSize',mrksz2,'DisplayName','Periodic extension')
    plot(xFp2D,yFp2D,'b.','MarkerSize',mrksz2,'HandleVisibility','off')
    
    xticks([0 pi/2 pi 3*pi/2 2*pi 2*pi+nxsupp/2*dxF])
    yticks([0 pi/2 pi 3*pi/2 2*pi 2*pi+nysupp/2*dyF])
    xticklabels({'0', '\pi/2','\pi','3\pi/2','2\pi','2\pi+w_x/2'})
    yticklabels({'0', '\pi/2','\pi','3\pi/2','2\pi','2\pi+w_y/2'})
    daspect([1 1 1])
    grid on