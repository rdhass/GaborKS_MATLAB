function plotGrid(xLES,yLES,zLES,nxQH,nyQH,nzQH,xQH,yQH,zQH,xF,yF,zF,xFp,yFp,gmxloc,gmyloc,gmzloc,nxsupp,nysupp)
    global gmxloc gmyloc gmzloc
    
    nxLES = length(xLES);
    nyLES = length(yLES);
    nzLES = length(zLES);
    xLES2D = xLES'*ones(1,nyLES);
    yLES2D = ones(nxLES,1)*yLES;

    xQH2D = xQH'*ones(1,nyQH+1);
    yQH2D = ones(nxQH+1,1)*yQH;
    
    nxF = length(xF);
    nyF = length(yF);
    nzF = length(zF);
    xF2D = xF'*ones(1,nyF);
    yF2D = ones(nxF,1)*yF;
    dxF = xF(2)-xF(1);
    dyF = yF(2)-yF(1);
    
    nxFp = length(xFp);
    nyFp = length(yFp);
    xFp2D = xFp'*ones(1,nyFp);
    yFp2D = ones(nxFp,1)*yFp;
    
    gmx = reshape(gmxloc,nxQH,nyQH,nzQH,[]);
    gmy = reshape(gmyloc,nxQH,nyQH,nzQH,[]);
    gmz = reshape(gmzloc,nxQH,nyQH,nzQH,[]);
     
    figure
    subplot(1,2,1), hold on
    plot(xLES2D(1,1),yLES2D(1,1),'rx','DisplayName','LES nodes')
    plot(xQH2D(1,1),yQH2D(1,1),'bo','DisplayName','QH boundary')
    plot(xF2D(1,1),yF2D(1,1),'g.','DisplayName','Fine grid')
    plot(squeeze(gmx(4,4,4,1)),squeeze(gmy(4,4,4,1)),'k.','DisplayName','Gabor mode locations')
    plot(xLES2D,yLES2D,'rx','HandleVisibility','off')
    plot(xQH2D,yQH2D,'bo','HandleVisibility','off')
    plot(xF2D,yF2D,'g.','HandleVisibility','off')
    plot(squeeze(gmx(4,4,4,:)),squeeze(gmy(4,4,4,:)),'k.','HandleVisibility','off')
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
    sgtitle('$x$-$y$ plane','Interpreter','Latex')
    
    xLES2D = xLES'*ones(1,nzLES);
    zLES2D = ones(nxLES,1)*zLES;
    xQH2D = xQH'*ones(1,nzQH+1);
    zQH2D = ones(nxQH+1,1)*zQH;
    xF2D = xF'*ones(1,nzF);
    zF2D = ones(nxF,1)*zF;
    
    figure, hold on
    plot(xLES2D(1,1),zLES2D(1,1),'rx','DisplayName','LES nodes')
    plot(xQH2D(1,1),zQH2D(1,1),'bo','DisplayName','QH boundary')
    plot(xF2D(1,1),zF2D(1,1),'g.','DisplayName','Fine grid')
    plot(squeeze(gmx(1,1,1,1)),squeeze(gmz(1,1,1,1)),'k.','DisplayName','Gabor mode locations')
    plot(xLES2D,zLES2D,'rx','HandleVisibility','off')
    plot(xQH2D,zQH2D,'bo','HandleVisibility','off')
    plot(xF2D,zF2D,'g.','HandleVisibility','off')
    plot(squeeze(gmx(1,1,1,:)),squeeze(gmz(1,1,1,:)),'k.','HandleVisibility','off')
    xticks([0 pi/2 pi 3*pi/2 2*pi])
    xticklabels({'0', '\pi/2','\pi','3\pi/2','2\pi'})
    grid on
    legend('Interpreter','Latex')
    daspect([1 1 1])
    
  
    title('$x$-$z$ plane','Interpreter','Latex')