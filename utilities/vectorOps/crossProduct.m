function [c1,c2,c3] = crossProduct(a1,a2,a3,b1,b2,b3)
    % Inputs:
    %  a1, a2, a3 --> Vector components. Each is a nxQH X nyQH X nzQH X nk X ntheta array
    %  b1, b2, b3 --> Vector components. Each is a nxQH X nyQH X nzQH X nk X ntheta array
    % Outputs:
    %  c1, c2, c3 --> Vector components such that c = a X b. Each is same
    %                 size as a and b vector components
    
    c1 = a2.*b3 - a3.*b2;
    c2 = a3.*b1 - a1.*b3;
    c3 = a1.*b2 - a2.*b1;