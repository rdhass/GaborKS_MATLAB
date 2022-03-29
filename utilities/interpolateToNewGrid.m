function [uout,vout,wout] = interpolateToNewGrid(uin,vin,win,xin,yin,zin,xout,yout,zout)
    if size(xin,1) == 1
        xin = transpose(xin);
    end
    if size(xout,1) == 1
        xout = transpose(xout);
    end
    uout = interp3(xin,yin,zin,uin,xout,yout,zout,'spline');
    vout = interp3(xin,yin,zin,vin,xout,yout,zout,'spline');
    wout = interp3(xin,yin,zin,win,xout,yout,zout,'spline');
    