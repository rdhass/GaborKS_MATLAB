function [norm1,norm2,norm3] = normalizeVec(a1,a2,a3)
    % Inputs:
    %  a1, a2, a3 --> Vector components. Each is a nxQH X nyQH X nzQH X nk X ntheta array
    % Outputs:
    %  norm1, norm2, norm3 --> Vector components such that norm = a/|a|
    
    amag = sqrt(a1.^2 + a2.^2 + a3.^2);
    norm1 = a1./amag;
    norm2 = a2./amag;
    norm3 = a3./amag;