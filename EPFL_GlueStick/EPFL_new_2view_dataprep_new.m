%%%%% EPFL_new two view estimations from GC_ransac
% addpath(genpath(pwd))
clear all;
 
rng('default')
s = rng;

% dataset = "FountainP11";
% dataset = "CastleP30";
dataset = "EntryP10";
% dataset = "HerzP8";
% dataset = "HerzP8";
dataset_name = "Entry_P10";

[K, R_true, t_true, im_names] = EPFL_get_ground_truth(dataset);
t_scale = sqrt(size(t_true,2))/norm(t_true);


%% Official Part of Data
% load("EPFL_code/EPFL_"+dataset+"_FCC_matched_full.mat")
estimated_fundamental_matrices = readmatrix(dataset_name+"_fundmats.txt");
n = size(estimated_fundamental_matrices,1)/3;
path_to_data=strcat('EPFL_code/Data/',dataset,'/',dataset,'_dense/');

% two_view_corresp = matched_pairs_epfl;
two_view_corresp = get_all_two_view_corresp(dataset_name, n);



% get ground truth emat and estimate emat
[Ks, Rs, ts, ccams, Tgt, imsizes] = EPFL_ground_truth_imsizes(im_names, path_to_data);
gtessenmats = calculate_gt_essenmats(ccams, n);
cam_intrinsics = calculate_intrinsics(Ks, imsizes, n);


essenmats = cell(n,n);
% inliers_idx = cell(n,n);
for i = 1:n-1
    for j = i+1:n
        if size(two_view_corresp{i,j},1) > 16
            cur_essen = Ks(:,:,2)' * estimated_fundamental_matrices(3*(i-1)+1:3*i,3*(j-1)+1:3*j) * Ks(:,:,1);
            cur_essen = cur_essen / norm(cur_essen);
            % [cur_essen, inliersIndex] = estimateEssentialMatrix(two_view_corresp{i,j}(:,[1,2]), two_view_corresp{i,j}(:,[3,4]), cam_intrinsics{i}, cam_intrinsics{j}, 'MaxNumTrials', 10000);
            % inliers_idx{i,j} = inliersIndex;

            essenmats{j,i} = cur_essen;
            essenmats{i,j} = essenmats{j,i}';
            
        end
    end
    fprintf("Done with estimating essential matrices for %d\n", i);
end

%%% Compare the estimated essential matrices and the ground truth 
% essenmat_metric = Compare_Essential_Matrices(gtessenmats, essenmats, n, two_view_corresp);

[E_est, data_keep] = emat_cell_to_matrix(essenmats);
[E, ~] = emat_cell_to_matrix(gtessenmats);

%%%% CHECKING SIGN FLIP HERE
% for i = 1:n
%     for j = 1:n
%         if size(essenmats{i,j},1) > 0
%             signflip = sign(sum(sign(essenmats{i,j} .* E(3*(i-1)+1:3*i,3*(j-1)+1:3*j)), "all"));
%             essenmats{i,j} = essenmats{i,j} * signflip;
%         end
%     end
% end
%%%% CHECKING SIGN FLIP HERE


% Prepare the data for LUD format
data.t = ts;
data.R = Rs;
data.n = n;
data.nPairs = 1;
data.keep = logical(data_keep);
data.E = E;
data.E_gt = E;
data.AdjMat = sparse(data_keep);
data.G_gt = Rs_reformatting(Rs);
data.E_est = essenmats;
[Hmat, tijGT] = get_estimated_relative_pose(essenmats, two_view_corresp, cam_intrinsics);


data.Hmat = Hmat;
data.tijGT = tijGT;
corr = cell(n,n);
for i = 1:n-1
    for j = i+1:n
        if size(two_view_corresp{i,j},1) > 10
            corr{i,j} = zeros(size(two_view_corresp{i,j},2)/2, size(two_view_corresp{i,j},1), 2);
            corr{i,j}(:,:,1) = two_view_corresp{i,j}(:,1:2)';
            corr{i,j}(:,:,2) = two_view_corresp{i,j}(:,3:4)';
        end
    end
end

data.corr = corr;
data.K = Ks;

Focal_gt = zeros(1,n);
imgW = zeros(n,1);
imgL = zeros(n,1);

for i = 1:n
    Focal_gt(i) = mean([Ks(1,1,i), Ks(2,2,i)]);
    imgW(i) = imsizes(1,1,i);
    imgL(i) = imsizes(2,1,i);
end
data.Focal_gt = Focal_gt;
data.imgW = imgW;
data.imgL = imgL;
data.W_mask_full = zeros(2,n);
data.W_x_full = zeros(2,n);
data.W_y_full = zeros(2,n);
data.CompInds = int32([1:n]');

data = gt_filter_pairs(data);


[rotation_metrics, angle_metrics] = Compare_estimated_gt_relative_rotations(data.G_gt, data.Hmat);
subplot(1,2,1)
histogram(rotation_metrics(rotation_metrics>0), 50)
title("rotation norm difference")
subplot(1,2,2)
histogram(angle_metrics(angle_metrics>0), 50)
title("rotation angle difference")

save("EPFL_GlueStick/"+dataset+"_two_view_data.mat", "data");




function [rotation_metrics, angle_metrics] = Compare_estimated_gt_relative_rotations(G_gt, Hmat)
n = size(Hmat,1)/3;
rotation_metrics = zeros(n,n);
angle_metrics = zeros(n,n);
for i = 1:n-1
    for j = i+1:n
        est_Rij = Hmat(3*(i-1)+1:3*i,3*(j-1)+1:3*j);
        if norm(est_Rij)>0
            gt_Rij = G_gt(3*(i-1)+1:3*i,3*(j-1)+1:3*j) / norm(G_gt(3*(i-1)+1:3*i,3*(j-1)+1:3*j));
            est_Rij = est_Rij / norm(est_Rij);
            rotation_metrics(i,j) = min(norm(est_Rij - gt_Rij), norm(est_Rij + gt_Rij));
            [angle_metrics(i,j),~] = new_CompareRotations(est_Rij, gt_Rij);
        end
    end
end
end

function [Hmat, tijGT] = get_estimated_relative_pose(essenmats, two_view_corresp, cam_intrinsics)

n = size(essenmats,1);
Hmat = zeros(3*n,3*n);
tijGT = cell(n,n);
for i = 1:n-1
    for j = i+1:n
        if size(essenmats{i,j},1)>0
            relativePose = estrelpose(essenmats{j,i}, cam_intrinsics{j}, cam_intrinsics{i}, two_view_corresp{i,j}(:,3:4), two_view_corresp{i,j}(:,1:2));
            Hmat(3*(i-1)+1:3*i,3*(j-1)+1:3*j) = relativePose.R;
            Hmat(3*(j-1)+1:3*j,3*(i-1)+1:3*i) = transpose(Hmat(3*(i-1)+1:3*i,3*(j-1)+1:3*j));
            tijGT{i,j} = -relativePose(1).Translation / norm(relativePose(1).Translation);
%             tijGT{j,i} = (-1)*relativePose.Translation;
        end
    end
end


end

function G_gt = Rs_reformatting(Rs)
n = size(Rs,3);
G_gt = zeros(3*n,3*n);
for i = 1:n
    for j = 1:n
        if i ~= j
            G_gt(3*(i-1)+1:3*i, 3*(j-1)+1:3*j) = Rs(:,:,i) * Rs(:,:,j)';
        end
    end
end

end

function [E_est, data_keep] = emat_cell_to_matrix(ematcell)
n = size(ematcell,1);
% E_est = cell(n,n);
E_est = zeros(3*n,3*n);
data_keep = zeros(n,n);
for i = 1:n-1
    for j = i+1:n
        if size(ematcell{j,i},1)>0 
            E_est(3*(j-1)+1:3*j, 3*(i-1)+1:3*i) = ematcell{j,i};
%             E_est{i,j} = ematcell{i,j};
%             E_est{j,i} = ematcell{i,j}';
            data_keep(i,j) = 1;
        end
    end
end

E_est = E_est + E_est';
data_keep = data_keep + data_keep';

end

function essenmat_metric = Compare_Essential_Matrices(gtessenmats, essenmats, n, two_view_corresp)

essenmat_metric = zeros(n,n);
for i = 1:n-1
    for j = i+1:n
        if size(two_view_corresp{i,j},1) > 10
            cur_gt_emat = gtessenmats{i,j}/norm(gtessenmats{i,j});
            cur_emat = essenmats{i,j}/norm(essenmats{i,j});
            essenmat_metric(i,j) = min(norm(cur_gt_emat - cur_emat), norm(cur_gt_emat + cur_emat));
        end        
    end
end
 
histogram(essenmat_metric(essenmat_metric >0),50)
title("Method Relative Essential_Matrix")

end

function gtessenmats = calculate_gt_essenmats(ccams,n)

gtessenmats = cell(n,n);
for i= 1:n-1
    for j = i+1:n
        cur_ess = vgg_F_from_P({ccams(:,:,i), ccams(:,:,j)});
        gtessenmats{j,i} = cur_ess / norm(cur_ess);
    end
end

end

function cam_intrinsics = calculate_intrinsics(Ks, imsizes, n)
    

cam_intrinsics = cell(1,n);

for i = 1:n
   cam_intrinsics{i} = cameraIntrinsics([Ks(1,1,i),Ks(2,2,i)], [Ks(1,3,i),Ks(2,3,i)], [imsizes(1,1,i), imsizes(2,1,i)]);  
end

end

function [elem] = get_distinct_elements(A, B)
% A and B are Nx4, Mx4 matrices, which correspond to
% keypoint locations. 
% Return: elem
%         elem: K x 4 matrix of A concatenated with the elements of B that
%         haven't appeared in A
if size(A, 1) == 0
    elem = B;
else
    [~, idx] = ismember(B, A, 'rows');
    bminusa = B(idx == 0, :);
    elem = [A; bminusa];
end

end






