
load("EPFL_code/"+dataset+"_sift_Z.mat");
EPFL_coords_data_final = SIFT;
n = length(SIFT);

S = FCC(Z, dimPerm', 2, 100, 8, 0);

matsize = dimPerm';


dimPerm_idx = [0, cumsum(matsize)]


matched_triplets_epfl = cell(n,n,n);
for i = 1:n
    for j = 1:n
        for k = 1:n
            if i < j && j < k
                matched_triplets_epfl{i,j,k} = [];
            end
        end
    end
end

% filter image correspondences
S = S > 0.99;

for i = 1:size(S,2)

    S_column = S(:,i);
    matchesnum = nnz(S_column);
   
    if matchesnum < 3
        continue
    end
    
    idx_data = zeros(2,matchesnum);
    idxs = find(S_column);
    for j = 1:length(idxs)
        cur_idx = sum(dimPerm_idx < idxs(j));

        feature_idx = idxs(j) - dimPerm_idx(cur_idx);

        idx_data(1,j) = cur_idx;
        idx_data(2,j) = feature_idx;
    end
    
    idx_data = idx_data';

    
    idx_data = sortrows(idx_data);

    d = size(idx_data,1);
    for j = 1:d-2
        for k = j+1:d-1
            for u = k+1:d
                cur_data = [EPFL_coords_data_final{idx_data(j,1)}.locs(idx_data(j,2),:), EPFL_coords_data_final{idx_data(k,1)}.locs(idx_data(k,2),:),EPFL_coords_data_final{idx_data(u,1)}.locs(idx_data(u,2),:)];
                    matched_triplets_epfl{idx_data(j,1),idx_data(k,1),idx_data(u,1)} = [matched_triplets_epfl{idx_data(j,1),idx_data(k,1),idx_data(u,1)}; cur_data];
            end
        end
    end
        
end
    
save("EPFL_code/EPFL_"+dataset+"_FCC_matched_triplets_full.mat", "matched_triplets_epfl", "S")


%% matching pairs
matched_pairs_epfl = cell(n,n);
for i = 1:n
    for j = 1:n
        if i < j
            matched_pairs_epfl{i,j} = [];
        end
    end
end

for i = 1:size(S,2)

    S_column = S(:,i);
    matchesnum = nnz(S_column);
   
    if matchesnum < 2
        continue
    end
    
    idx_data = zeros(2,matchesnum);
    idxs = find(S_column);
    for j = 1:length(idxs)
        cur_idx = sum(dimPerm_idx < idxs(j));

        feature_idx = idxs(j) - dimPerm_idx(cur_idx);

        idx_data(1,j) = cur_idx;
        idx_data(2,j) = feature_idx;
    end
    
    idx_data = idx_data';
         

    idx_data = sortrows(idx_data);

    d = size(idx_data,1);
    for j = 1:d-1
        for k = j+1:d
            cur_data = [EPFL_coords_data_final{idx_data(j,1)}.locs(idx_data(j,2),:), EPFL_coords_data_final{idx_data(k,1)}.locs(idx_data(k,2),:)];
            matched_pairs_epfl{idx_data(j,1),idx_data(k,1)} = [matched_pairs_epfl{idx_data(j,1),idx_data(k,1)}; cur_data];
        end
    end
        
end

save("EPFL_code/EPFL_"+dataset+"_FCC_matched_full.mat")






