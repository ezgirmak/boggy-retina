function [info] = CompartmentInfo(tree, index)
    
    compType = typeN_tree(tree);
    idpar = idpar_tree(tree);
    seglength = len_tree(tree);
    
    % Averages for cell for comparison
    dx_avg = mean(abs(tree.X - tree.X(idpar)));
    dy_avg = mean(abs(tree.Y - tree.Y(idpar)));
    dz_avg = mean(abs(tree.Z - tree.Z(idpar)));
    D_avg = mean(tree.D);
    seglen_avg = mean(seglength);
    Ri_avg = mean(tree.Ri);

    indices = [index; find(tree.dA(index, :)); find(tree.dA(:, index))];
    dx = [dx_avg; abs(tree.X(index) - tree.X(indices))];
    dy = [dy_avg; abs(tree.Y(index) - tree.Y(indices))];
    dz = [dz_avg; abs(tree.Z(index) - tree.Z(indices))];
    D = [D_avg; tree.D(indices)];
    seglen = [seglen_avg; seglength(indices)];
    Ri = [Ri_avg; tree.Ri(indices)];
    R = [-1; tree.R(indices)];
    type = [-1; compType(indices)];
    indices = [0; indices];
    
    info = table(indices, dx, dy, dz, D, seglen, ...
        R, Ri, type);

end