function fprime = ddzFD(f,dz)
    fprime = zeros(size(f));
    fprime(:,:,2:end-1) = (f(:,:,3:end) - f(:,:,1:end-2))./(2*dz);
    fprime (:,:,1) = (-3*f(:,:,1) + 4*f(:,:,2) - f(:,:,3))./(2*dz);
    fprime(:,:,end) = (3*f(:,:,end) - 4*f(:,:,end-1) + f(:,:,end-2))./(2*dz);