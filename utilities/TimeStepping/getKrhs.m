function rhs = getKrhs(k, dudx)
    rhs = zeros(3,1);
    rhs(1) = -(k(1)*dudx(1,1) + k(2)*dudx(2,1) + k(3)*dudx(3,1));
    rhs(2) = -(k(1)*dudx(1,2) + k(2)*dudx(2,2) + k(3)*dudx(3,2));
    rhs(3) = -(k(1)*dudx(1,3) + k(2)*dudx(2,3) + k(3)*dudx(3,3));