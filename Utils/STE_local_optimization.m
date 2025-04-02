% inputs should be T, Corresp N x 6, CalM, and a threshold (if threshold 
% is too small, at least 8 points will be used

function [T, distances2, inlier_Corresp, R_t_2, R_t_3] = STE_local_optimization(T, Corresp, distances, thresh_percentage, CalM)

num_matches = size(distances,2);
subspace_dist_means = mean(reshape(distances,[9, num_matches/9]),1);
% threshold change to percentage, so assume that the number of inliers is
% around 50 percent

num_inlier_thresh = round(size(subspace_dist_means,2) * thresh_percentage);
[new_distances,inlier_indices] = mink(subspace_dist_means, min(max(num_inlier_thresh, 8),30));

% % inlier_indices are the indices of the smallest threshold distances / inliers. Inlier_Corresp are the inlier data
inlier_Corresp = Corresp(inlier_indices, :);

% % Using the inliers, we estimate the trifocal tensor using more than 7 data, in hope that it gives a better estimate. 
[R_t_2,R_t_3,Reconst,cur_Test2,iter] = ResslTFTPoseEstimation(inlier_Corresp', CalM);
% [R_t_2,R_t_3,Reconst,cur_Test2,iter] = LinearTFTPoseEstimation(inlier_Corresp', CalM);

T = cur_Test2;

distances2 = ransac_ReprError({CalM(1:3,:) * eye(3,4),CalM(4:6,:) *R_t_2, CalM(7:9,:) * R_t_3}, inlier_Corresp');

end

