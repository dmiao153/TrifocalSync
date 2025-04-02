function np = reorder_np_array_to_cell(x)

    n = size(x,1)/3;
    np = {};
    
    for i = 1:n 
        np{i} = x(3*(i-1)+1:3*i,:);
    end
    
    
end

