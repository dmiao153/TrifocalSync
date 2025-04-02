function data = gt_filter_pairs(data)
% This function filters out the pairs of images using ground truth, where
% the orientation of their view is completely different. 

n = data.n;

rotation_differences = zeros(n,n);
for i = 1:n
    for j = 1:n
        [E,e,normE,norme]=CompareRotations(data.R(:,:,i),data.R(:,:,j));
        rotation_differences(i,j) = e;
    end
end

new_keep = data.keep;
for i = 1:n
    for j = 1:n
        if rotation_differences(i,j) > 80
            new_keep(i,j) = 0;
            data.E(3*(i-1)+1:3*i, 3*(j-1)+1:3*j) = zeros(3,3);
            data.E_gt(3*(i-1)+1:3*i, 3*(j-1)+1:3*j) = zeros(3,3);
            data.E_est{i,j} = zeros(3,3);
            data.Hmat(3*(i-1)+1:3*i, 3*(j-1)+1:3*j) = zeros(3,3);
            data.G_gt(3*(i-1)+1:3*i, 3*(j-1)+1:3*j) = zeros(3,3);

            data.tijGT{i,j} = [];
        end
    end
end

data.AdjMat = sparse(new_keep);
end