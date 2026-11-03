function G = dipole_bz_kernel(X, Y, position)
%DIPOLE_BZ_KERNEL Bz in nT per magnetic moment in A*m^2.
% position = [x_m, y_m, depth_m], with depth POSITIVE BELOW the z=0 plane.
% Magnetization components use a right-handed frame with z POSITIVE UP.
% Bz(:) = G * [mx;my;mz]. This is a point-dipole forward model.
    assert(isequal(size(X),size(Y)), 'X and Y must have the same size.');
    assert(numel(position)==3 && position(3)>0, 'Depth must be positive.');
    dx = X(:)-position(1);
    dy = Y(:)-position(2);
    dz = position(3);
    r2 = dx.^2+dy.^2+dz.^2;
    r5 = r2.^2.5;
    G = 100*[3*dz*dx./r5, 3*dz*dy./r5, (3*dz^2-r2)./r5];
end
