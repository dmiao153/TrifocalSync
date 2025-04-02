function [final_data, Filter_ind] = prune_sparse_camera_blocks(final_data)
n = size(final_data,1);
w_mask = zeros(n,n,n);

for i = 1:n
    for j = 1:n
        for k = 1:n 
            if i~=j || j~=k || i~=k 
                if size(final_data{i,j,k},1) > 0
                    w_mask(i,j,k) = 1;
                end
            else
                w_mask(i,i,i) = 1;
            end         
        end
    end
end


% If we observed more than 5 percent of all the entries, then we keep the
% camera. Otherwise, we will prune it once again. 
Filter_ind = [];
for i = 1:n
    if nnz(w_mask(:,:,i)) > 0.05 * n^2
        Filter_ind = [Filter_ind, i];
    end
end

final_data = final_data(Filter_ind, Filter_ind, Filter_ind);

end