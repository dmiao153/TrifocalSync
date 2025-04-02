function [indices, cm, Ps, check, uPs, trifocal_tensors] = getCameras(T)
    n = size(T,1)/3;
    
    trifocal_tensors = {};
    indices = zeros(n-3, 6);
    
    for i = 1: (n-3)
        trifocal_tensors{i,1} = {retrieve_individual_trifocal_tensor(T,i,i+1,i+2)};
        trifocal_tensors{i,2} = {retrieve_individual_trifocal_tensor(T,i+1,i+2,i+3)};
        indices(i,:) = [i,i+1,i+2,i+1,i+2,i+3];
    end
        
    [indices, cm, Ps, check, uPs] = getCameras_from_block_tft(indices, trifocal_tensors, n);

end

