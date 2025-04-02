function npmat = reorder_np_cell_to_array(np)

    n = max(size(np));
    npmat = zeros(3*n,4);
    for i = 1:n
        npmat(3*(i-1)+1:3*i, :) = np{i};
    end
    
end