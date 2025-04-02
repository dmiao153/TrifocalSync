%%% This script calculates trifocal tensors linearly from line information
%%% and double checks it against the ground truth information. 
rng('default')
s = rng;

dataset = "HerzP8";
% dataset = "FountainP11";
dataset_name = "Herz_P8";
n = 8;

load("EPFL_GlueStick/"+dataset+"_two_view_data.mat");

n = data.n; 

for i = 1:n
    for j = 1:n
        if nnz(data.E_est{i,j}) > 0
            data.E(3*(i-1)+1:3*i, 3*(j-1)+1:3*j) = data.E(3*(i-1)+1:3*i, 3*(j-1)+1:3*j) / norm(data.E(3*(i-1)+1:3*i, 3*(j-1)+1:3*j));
            data.E_est{i,j} = data.E_est{i,j} / norm(data.E_est{i,j});
        else
            data.E(3*(i-1)+1:3*i, 3*(j-1)+1:3*j) = zeros(3,3);
            data.corr{i,j} = [];
        end
    end
end
data.E_gt = data.E;



[triplet_points, triplet_lines, n] = EPFL_process_data(dataset_name, n);


[K, R_true, t_true, im_names] = EPFL_get_ground_truth(dataset);
t_scale = sqrt(size(t_true,2))/norm(t_true);


path_to_data=strcat('EPFL_code/Data/',dataset,'/',dataset,'_dense/');
dir_data = strcat('EPFL_code/Data/',dataset,'/images/');

initial_sample_size=100;
bundle_adj_size=50;
repr_err_th=1;

[Ks, Rs, ts, ccams, cTgt] = EPFL_ground_truth(im_names, path_to_data, t_scale);

% load("EPFL_code/EPFL_"+dataset+"_FCC_matched_full.mat")
% corresp_by_triplet = triplet_points;
% n = size(corresp_by_triplet,1);

% tft_estimated_cams = cell(n,n,n);
% data_count = zeros(n,n,n);
cTest = cell(n,n,n);
for i= 1:n
    for j = 1:n
        for k = 1:n
            if i < j && j < k && nnz(data.E(3*(j-1)+1:3*j,3*(i-1)+1:3*i)) > 0 ...
                && nnz(data.E(3*(k-1)+1:3*k,3*(j-1)+1:3*j)) > 0 ...
                && nnz(data.E(3*(k-1)+1:3*k,3*(i-1)+1:3*i)) > 0 ...
                && size(triplet_points{i,j,k},1) > 10
                
                % [P1,P2,P3] = tft_from_3_essmat(data.E_est{i,j},...
                %     data.E_est{i,k}, ...
                %     data.E_est{j,k});
                
%                 [P1,P2,P3] = tft_from_3_essmat(data.E_est{j,i},...
%                                     data.E_est{k,i}, ...
%                                     data.E_est{k,j});
                
%                 [P1,P2,P3] = tft_from_3_essmat(data.E(3*(j-1)+1:3*j,3*(i-1)+1:3*i),...
%                     data.E(3*(k-1)+1:3*k,3*(i-1)+1:3*i), ...
%                     data.E(3*(k-1)+1:3*k,3*(j-1)+1:3*j));
                
%                 [P1,P2,P3] = tft_from_3_essmat(data.E(3*(j-1)+1:3*j,3*(i-1)+1:3*i)',...
%                     data.E(3*(k-1)+1:3*k,3*(i-1)+1:3*i)', ...
%                     data.E(3*(k-1)+1:3*k,3*(j-1)+1:3*j)');
                
%                             if i < j && j < k && nnz(data.E_est{j,i}) > 0 ...
%                 && nnz(data.E_est{k,j}) > 0 ...
%                 && nnz(data.E_est{k,i}) > 0
%                 
                [P1,P2,P3] = tft_from_3_essmat(data.E_est{j,i},...
                    data.E_est{k,i}, ...
                    data.E_est{k,j});



                tempt = T_from_P({P1,P2,P3});
                [R_t_2, R_t_3] = R_t_from_TFT(tempt, [eye(3);eye(3);eye(3)], triplet_points{i,j,k}');
                 
                N = size(triplet_points{i,j,k},1);
 
                init_sample=randsample(1:N,min(50,N));
 
                [R_t_ref,~,nit,repr_errBA]=BundleAdjustment([K;K;K],...
                    [eye(3,4);R_t_2;R_t_3],triplet_points{i,j,k}(init_sample,:)');

                 
                % ref_R_t = {R_t_ref(1:3,:),R_t_ref(4:6,:), R_t_ref(7:9,:)};
                im1 = i; im2= j; im3= k;
                temp231 = T_from_P({R_t_ref(1:3,:),R_t_ref(4:6,:), R_t_ref(7:9,:)});
                temp321 = T_from_P({R_t_ref(1:3,:),R_t_ref(7:9,:), R_t_ref(4:6,:)});
                temp132 = T_from_P({R_t_ref(4:6,:),R_t_ref(1:3,:), R_t_ref(7:9,:)});
                temp312 = T_from_P({R_t_ref(4:6,:),R_t_ref(7:9,:), R_t_ref(1:3,:)});
                temp123 = T_from_P({R_t_ref(7:9,:),R_t_ref(1:3,:), R_t_ref(4:6,:)});
                temp213 = T_from_P({R_t_ref(7:9,:),R_t_ref(4:6,:), R_t_ref(1:3,:)});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                 tempt = T_from_P({P1,P2,P3});
%                 [a1,a2,a3] = P_from_T(tempt/frob(tempt));
%                 [a2,a3] = calibrate(a2,a3);
%                 cTest{j,k,i} = T_from_P({a1,a2,a3});
%                 
%                 
%                 tempt = T_from_P({P1,P3,P2});
%                 [a1,a2,a3] = P_from_T(tempt/frob(tempt));
%                 [a2,a3] = calibrate(a2,a3);
%                 cTest{k,j,i} = T_from_P({a1,a2,a3});
%                 
%                 
%                 tempt = T_from_P({P2,P1,P3});
%                 [a1,a2,a3] = P_from_T(tempt/frob(tempt));
%                 [a2,a3] = calibrate(a2,a3);
%                 cTest{i,k,j} = T_from_P({a1,a2,a3});
%                 
%                 
%                 tempt = T_from_P({P2,P3,P1});
%                 [a1,a2,a3] = P_from_T(tempt/frob(tempt));
%                 [a2,a3] = calibrate(a2,a3);
%                 cTest{k,i,j} = T_from_P({a1,a2,a3});
%                 
%                 
%                 tempt = T_from_P({P3,P1,P2});
%                 [a1,a2,a3] = P_from_T(tempt/frob(tempt));
%                 [a2,a3] = calibrate(a2,a3);
%                 cTest{i,j,k} = T_from_P({a1,a2,a3});
%                 
%                 
%                 tempt = T_from_P({P3,P2,P1});
%                 [a1,a2,a3] = P_from_T(tempt/frob(tempt));
%                 [a2,a3] = calibrate(a2,a3);
%                 cTest{j,i,k} = T_from_P({a1,a2,a3});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                tempt = T_from_P({P1,P2,P3});
                gttt = temp231;
                % gttt = T_from_P({ccams(:,:,i), ccams(:,:,j), ccams(:,:,k)});
                sFlip = sign(sum(sign(tempt.*gttt),'all'));
                if sFlip == -1
                    disp("pause")
                end
                % tempt = transform_TFT(tempt, K,K,K, 0);
                % sFlip = fix_cheirality_new(tempt, triplet_points{i,j,k}');
%                 cdata = transform_to_calibrated(triplet_points{i,j,k}, K);
%                 sFlip = fix_cheirality_new(tempt, cdata');
                
%                 sFlip = 1;
                cTest{j,k,i} = sFlip * T_from_P({P1,P2,P3});
                % 
                % 
                tempt = T_from_P({P1,P3,P2});
                gttt = temp321;
                % gttt = T_from_P({ccams(:,:,i), ccams(:,:,k), ccams(:,:,j)});
                sFlip = sign(sum(sign(tempt.*gttt),'all'))  ;
                if sFlip == -1
                    disp("pause")
                end
                % tempt = transform_TFT(tempt, K,K,K, 0);
                % % % sFlip = fix_cheirality_new(tempt, triplet_points{i,j,k}(:,[1,2,5,6,3,4])');
                % % cdata = transform_to_calibrated(triplet_points{i,j,k}(:,[1,2,5,6,3,4]), K);
                % % sFlip = fix_cheirality_new(tempt, cdata');
                % % 
                cTest{k,j,i} = sFlip * T_from_P({P1,P3,P2});
               
                tempt = T_from_P({P2,P1,P3});
                gttt = temp132;
                % gttt = T_from_P({ccams(:,:,j), ccams(:,:,i), ccams(:,:,k)});
                sFlip = sign(sum(sign(tempt.*gttt),'all'));                % tempt = transform_TFT(tempt, K,K,K, 0);
                if sFlip == -1
                    disp("pause")
                end
                % % % sFlip = fix_cheirality_new(tempt, triplet_points{i,j,k}(:,[3,4,1,2,5,6])');
                % % cdata = transform_to_calibrated(triplet_points{i,j,k}(:,[3,4,1,2,5,6]), K);
                % % sFlip = fix_cheirality_new(tempt, cdata');
                % % 
                cTest{i,k,j} = sFlip * T_from_P({P2,P1,P3});
                % % 
%                 
                tempt = T_from_P({P2,P3,P1});
                gttt = temp312;
                % gttt = T_from_P({ccams(:,:,j), ccams(:,:,k), ccams(:,:,i)});
                sFlip = sign(sum(sign(tempt.*gttt),'all'));
                if sFlip == -1
                    disp("pause")
                end
                % % tempt = transform_TFT(tempt, K,K,K, 0);
                % % sFlip = fix_cheirality_new(tempt, triplet_points{i,j,k}(:,[3,4,5,6,1,2])');
                % % cdata = transform_to_calibrated(triplet_points{i,j,k}(:,[3,4,5,6,1,2]), K);
                % % sFlip = fix_cheirality_new(tempt, cdata');
 
                cTest{k,i,j} = sFlip * T_from_P({P2,P3,P1});
                % % 
                tempt = T_from_P({P3,P1,P2});
                % gttt = T_from_P({ccams(:,:,k), ccams(:,:,i), ccams(:,:,j)});
                gttt = temp123;
                sFlip = sign(sum(sign(tempt.*gttt),'all'));
                % % % tempt = transform_TFT(tempt, K,K,K, 0);
                % % % sFlip = fix_cheirality_new(tempt, triplet_points{i,j,k}(:,[5,6,1,2,3,4])');
                % % cdata = transform_to_calibrated(triplet_points{i,j,k}(:,[5,6,1,2,3,4]), K);
                % % sFlip = fix_cheirality_new(tempt, cdata');
                % % 
                if sFlip == -1
                    disp("pause")
                end
                cTest{i,j,k} = sFlip * T_from_P({P3,P1,P2});
                % % 
                % % 
                tempt = T_from_P({P3,P2,P1});
                % gttt = T_from_P({ccams(:,:,k), ccams(:,:,j), ccams(:,:,i)});
                gttt = temp213;
                sFlip = sign(sum(sign(tempt.*gttt),'all'));
%                 % % 
                % % % tempt = transform_TFT(tempt, K,K,K, 0);
                % % % sFlip = fix_cheirality_new(tempt, triplet_points{i,j,k}(:,[5,6,3,4,1,2])');
                % % cdata = transform_to_calibrated(triplet_points{i,j,k}(:,[5,6,3,4,1,2]), K);
                % % sFlip = fix_cheirality_new(tempt, cdata');
                % % 
                if sFlip == -1
                    disp("pause")
                end
                cTest{j,i,k} = sFlip * T_from_P({P3,P2,P1});

            end
        end
    end
    fprintf("iteration %d \n", i);
end

% for i = 1:n
%     for j = 1:n
%         if i ~= j
%             cTest{i,i,j} = reorder_iij_from_fundamental(data.E_est{i,j});
%         end
%     end
% end

[Ks, Rs, ts, ccams, cTgt] = EPFL_ground_truth(im_names, path_to_data, t_scale);



cameras_gt = cell(1,n);
for i = 1:n
   cameras_gt{i} = ccams(:,:,i); 
end


[metrics, metrics_rotation, metrics_translation] = compare_block_tensor(cTest, n, cameras_gt);
figure;
subplot(1,3,1)
histogram(metrics(metrics >0),50)
title("Method Relative TFT")

subplot(1,3,2)
metrics_rotation = metrics_rotation(metrics_rotation > 0);
histogram(metrics_rotation(:),50)
title("Method Rotation Errors")

subplot(1,3,3)
metrics_translation = metrics_translation(metrics_translation > 0);
histogram(metrics_translation(:),50)
title("Method Translation Errors")


T = final_T_cell_to_array(cTest);

% % %%% HARD CODE FOR HERZP25, AS IT IS THE ONLY ONE THAT HAS CAMERAS WITH
% % %%% MISSING MEASUREMENTS
% % % % % % % good_indices = [1,2,3,4,5,6,7,8,9,10,11,12,13,15,16,17,18,19,20,21,22,23,24,25];
% % % % % % % cTest2 = cTest(good_indices, good_indices, good_indices);
% % % % % % % cameras_gt2 = cameras_gt(good_indices);
% % % % % % % fprintf("We actually have only %d cameras out of 25, hardcoded to take out\n", length(good_indices));
% % % % % % % T2 = final_T_cell_to_array(cTest2);
% % % % % % % clear cTest2
% % % % % % % final_T = heuristic_hosvd_imputation(T2, cameras_gt2, t_scale);
tic
final_T = heuristic_hosvd_imputation_final(T, cameras_gt, t_scale, true);
toc
% save("EPFL_code/"+dataset+"_final_tensor.mat", "T", "cameras_gt", "final_T")


% load("EPFL_code/"+dataset+"_final_tensor.mat")


res = residual_hosvd(final_T, cameras_gt, t_scale);
disp('HO-rSTE');
disp('Relative Cam Difference | Mean Location Error | Median Location Error | Mean Rotation Error | Median Rotation Error |');
fprintf('%.6f %.6f %.6f %.6f %.6f \n', res.res2, res.errT_mean, res.errT_median, res.errR_mean, res.errR_median);

% save(dataset+"_cameras_est.mat")

function cdata = transform_to_calibrated(triplet_points, K)
n = size(triplet_points,1);

data1 = [triplet_points(:,1:2), ones(n,1)];
data2 = [triplet_points(:,3:4), ones(n,1)];
data3 = [triplet_points(:,5:6), ones(n,1)];
  
cdata1 = inv(K) * data1';
cdata2 = inv(K) * data2';
cdata3 = inv(K) * data3';
    
cdata1 = cdata1 ./ repmat(cdata1(3,:), [3,1]);
cdata2 = cdata2 ./ repmat(cdata2(3,:), [3,1]);
cdata3 = cdata3 ./ repmat(cdata3(3,:), [3,1]);

cdata1 = cdata1';
cdata2 = cdata2';
cdata3 = cdata3';

cdata = [cdata1(:,1:2), cdata1(:,1:2), cdata1(:,1:2)];

end



function res = residual_hosvd(cur_t, np, t_scale)

cameras = nvecs(sptensor(cur_t), 1, 4);

[dd2,nn2,Hf, npmat1, npmat2] = compare_projection_matrices(reorder_np_array_to_cell(cameras), np, false);
res2 = dd2 / nn2;

cameras = normalizenpmat4(cameras);
np = normalizenpmat4(reorder_np_cell_to_array(np));

[errR_mean, errR_median, errT_mean, errT_median] = comparison_for_calibrated_cameras(reorder_np_array_to_cell(cameras), reorder_np_array_to_cell(np));
res.errR_mean = errR_mean;
res.errR_median = errR_median;
res.errT_mean = errT_mean / t_scale;
res.errT_median = errT_median / t_scale;
res.res2 = res2;
end

function npmat = normalizenpmat4(npmat)
    n = size(npmat,1)/3;
%     for i = 1:n
%        npmat(3*(i-1)+1:3*i,4) = npmat(3*(i-1)+1:3*i,4) / norm(npmat(3*(i-1)+1:3*i,4),'fro'); 
%     end
    for i = 1:n
       npmat(3*(i-1)+1:3*i,:) = npmat(3*(i-1)+1:3*i,:) / norm(npmat(3*(i-1)+1:3*i,1:3),'fro'); 
    end
end

function [metrics, metrics_rotation, metrics_translation] = compare_block_tensor(cTest, n, cameras_gt)

T_gt = generate_block_trifocal_tensor(cameras_gt, n);

metrics = zeros(n,n,n);
metrics_rotation = zeros(n,n,n);
metrics_translation = zeros(n,n,n);

for i= 1:n
    for j = 1:n
        for k = 1:n
            if nnz(cTest{i,j,k}) > 0
                metrics(i,j,k) = scaled_diff_frobenius_norm(cTest{i,j,k}, T_gt(3*(i-1)+1:3*i, 3*(j-1)+1:3*j, 3*(k-1)+1:3*k));
                [errR_mean, errT_mean] = angle_diff(cTest{i,j,k}, T_gt(3*(i-1)+1:3*i, 3*(j-1)+1:3*j, 3*(k-1)+1:3*k));
                metrics_rotation(i,j,k) = errR_mean;
                metrics_translation(i,j,k) = errT_mean;
            end
        end
    end
end


end





% % % estimate_tft_for_bundler
% % % addpath(genpath(pwd))
% % dataset = "Piccadilly";
% % load("External Packages and Toolboxes/SfM_CVPR2017_code/Data/"+dataset+"_data.mat")
% % Jmat = diag([1 1 -1]);
% % n = data.n;
% % Rmat_orig = zeros(3,3,n);
% % tmat_orig = zeros(3,1,n);
% % ccams = zeros(3,4,n);
% % cameras_gt = cell(1,n);
% % t_scale = sqrt(n) / frob(data.t);
% % for i = 1:n
% %     Rmat_orig(:,:,i) = data.R(:,:,i)';
% %     tmat_orig(:,:,i) = - Rmat_orig(:,:,i) * data.t(:,i) * t_scale;
% %     ccams(:,:,i) = Rmat_orig(:,:,i)' * [eye(3), -tmat_orig(:,:,i)];
% %     cameras_gt{i} = ccams(:,:,i);
% % end
% % 
% % good_index = [];
% % for i = 1:n
% %     Rmat_orig(:,:,i) = data.R(:,:,i)';
% %     tmat_orig(:,:,i) = - Rmat_orig(:,:,i) * data.t(:,i) * t_scale;
% %     ccams(:,:,i) = Rmat_orig(:,:,i)' * [eye(3), -tmat_orig(:,:,i)];
% %     cameras_gt{i} = ccams(:,:,i);
% %     if nnz(cameras_gt{i}) > 0
% %         good_index = [good_index, i];
% %     end
% % end
% % randsam = randsample(length(good_index), 1000);
% % good_index = good_index(randsam);
% % load("Photo_Tourism_cleaning/" + dataset+"/"+dataset+"good_index.mat")
% % % save("Photo_Tourism_cleaning/" + dataset+"/"+dataset+"good_index.mat", "good_index")
% % 
% % %% filter the void estimates
% % E_est = data.E_est(good_index, good_index);
% % cameras_gt = cameras_gt(good_index);
% % ccams = ccams(:,:,good_index);
% % Rmat_orig = Rmat_orig(:,:,good_index);
% % tmat_orig = tmat_orig(:,:,good_index);
% % n = length(cameras_gt);
% % 
% % fprintf("Just cleared all the void initial estimates\n")
% % clear data;
% % 
% % %% calculate cameras from essential matrices
% % tft_estimated_cams = cell(n,n,n);
% % data_count = zeros(n,n,n);
% % parfor i= 1:n
% %     for j = 1:n
% %         for k = 1:n
% %             if i < j && j < k && size(E_est{j,i},1) > 0 && size(E_est{k,j},1) > 0 && size(E_est{k,i},1) > 0
% % 
% %                 [P1,P2,P3] = tft_from_3_essmat(E_est{j,i}, E_est{k,i}, E_est{k,j});
% %                 tft_estimated_cams{i,j,k} = {P1,P2,P3};
% %                 data_count(i,j,k) = 1;
% %             end
% %         end
% %     end
% %     fprintf("iteration %d \n", i);
% % end
% % 
% % for i = 1:n
% %     for j = 1:n
% %         for k = 1:n
% %             if i<j && j < k && data_count(i,j,k)
% %                 data_count(j,j,k) = 1;
% %                 data_count(k,k,j) = 1;
% %                 data_count(i,i,j) = 1;
% %                 data_count(j,j,i) = 1;
% %                 data_count(i,i,k) = 1;
% %                 data_count(k,k,i) = 1;
% % 
% %                 data_count(j,k,i) = 1;
% %                 data_count(k,j,i) = 1;
% %                 data_count(i,k,j) = 1;
% %                 data_count(j,i,j) = 1;
% %                 data_count(j,i,k) = 1;
% %             end
% %         end
% %     end
% % end
% % 
% % %% We first truncate the size of cTest then calculate it 
% % datacounts = sum(sum(data_count));
% % first_filtering = datacounts > (0.015 * n^2);
% % first_filtering = find(first_filtering);
% % if length(first_filtering) > 225
% %     [~,first_filtering]  = maxk(datacounts, 225);
% % end
% % save("Photo_Tourism_cleaning/" + dataset+"/"+dataset+"first_filtering_indices.mat", "first_filtering")
% % 
% % fprintf("The size of the first filtering is %d\n", length(first_filtering));
% % tft_estimated_cams = tft_estimated_cams(first_filtering, first_filtering, first_filtering);
% % cameras_gt = cameras_gt(first_filtering);
% % E_est = E_est(first_filtering, first_filtering);
% % n = length(first_filtering);
% % data_count = data_count(first_filtering, first_filtering, first_filtering);
% 
% load("Photo_Tourism_cleaning/Piccadilly/Piccadilly_newindices.mat")
% fprintf("The size of newinds is %d\n", length(newinds));
% tft_estimated_cams = tft_estimated_cams(newinds, newinds, newinds);
% cameras_gt = cameras_gt(newinds);
% E_est = E_est(newinds, newinds);
% n = length(newinds);
% 
% 
% cTest = cell(n,n,n);
% for i = 1:n
%     for j = 1:n
%         for k = 1:n
%             if i <j && j < k && size(tft_estimated_cams{i,j,k}, 1)>0
%                 P1 = tft_estimated_cams{i,j,k}{1};
%                 P2 = tft_estimated_cams{i,j,k}{2};
%                 P3 = tft_estimated_cams{i,j,k}{3};
% 
% 
%                 cTest{j,k,i} = T_from_P({P1,P2,P3});
%                 cTest{k,j,i} = T_from_P({P1,P3,P2});
%                 cTest{i,k,j} = T_from_P({P2,P1,P3});
%                 cTest{k,i,j} = T_from_P({P2,P3,P1});
%                 cTest{i,j,k} = T_from_P({P3,P1,P2});
%                 cTest{j,i,k} = T_from_P({P3,P2,P1});
% 
%             end
%         end
%     end
% end
% clear tft_estimated_cams datacounts data_count;
% % save("Photo_Tourism_cleaning/"+dataset+"/"+dataset+"_tft_from_essmat.mat", "cTest", '-v7.3');
% 
% fprintf("loading the estimated trifocal tensors \n");
% % load("Photo_Tourism_cleaning/"+dataset+"/"+dataset+"_tft_from_essmat.mat")
% parfor i = 1:n
%     for j = 1:n
%         if nnz(E_est{i,j}) > 0
%             cTest{j,j,i} = reorder_iij_from_fundamental(E_est{i,j}); 
%         end
%     end
% end
% 
% 
% diary Photo_Tourism_cleaning/Piccadilly/output.log
% fprintf("Total number of cameras is %d \n", length(first_filtering));
% try
%     fprintf("%s starting to compare trifocal tensors \n", dataset);
% catch
% end
% 
% 
% T = final_T_cell_to_array(cTest);
% clear cTest;
% fprintf("Starting to synchronize\n");
% tic
% final_T = heuristic_hosvd_imputation(T, cameras_gt, t_scale);
% toc
% save("Photo_Tourism_cleaning/"+dataset+"/"+dataset+"_final_tensor_filtered.mat","cameras_gt", "final_T", '-v7.3')
% 
% 
% % load("EPFL_code/"+dataset+"_final_tensor.mat")
% 
% [res, cameras_est] = residual_hosvd(final_T, cameras_gt, t_scale);
% disp('HOSVD');
% disp('Relative Cam Difference | Mean Location Error | Median Location Error | Mean Rotation Error | Median Rotation Error |');
% fprintf('%.4f %.4f %.4f %.4f %.4f \n', res.res2, res.errT_mean, res.errT_median, res.errR_mean, res.errR_median);
% diary off
% 
% save("Photo_Tourism_cleaning/"+dataset+"/"+dataset+"_cameras_est_filtered.mat", "res", "cameras_est", '-v7.3');
% 
% 
% function cams_good = get_cams_from_final_T(final_T, good_index)
% [U,S,sv] = mlsvd(final_T);
% cs = U{1}(:,1:4);
% B = U{2}(:,1:4);
% C = U{3}(:,1:6);
% G = S(1:4, 1:4, 1:6);
% cur_t2 = tmprod(G, {cs,B,C}, 1:3);
% 
% 
% cur_t2(isnan(cur_t2)) = 0;
% cameras = retrieve_cameras_hosvd(cur_t2);
% 
% cams_good = zeros(3*length(good_index), 4);
% for i = 1:length(good_index)
%    cams_good(3*(i-1)+1:3*i,:) = cameras(3*(i-1)+1:3*i, :); 
% end
% 
% end
% 
% function [res, cameras_est] = residual_hosvd(cur_t, np, t_scale)
% [U,S,sv] = mlsvd(cur_t);
% cs = U{1}(:,1:4);
% B = U{2}(:,1:4);
% C = U{3}(:,1:6);
% G = S(1:4, 1:4, 1:6);
% cur_t2 = tmprod(G, {cs,B,C}, 1:3);
% 
% 
% cur_t2(isnan(cur_t2)) = 0;
% cameras = retrieve_cameras_hosvd(cur_t2);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uncali = false;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [dd2,nn2,Hf, npmat1, npmat2] = compare_projection_matrices(reorder_np_array_to_cell(cameras), np, uncali);
% res2 = dd2 / nn2;
% 
% cameras = normalizenpmat4(cameras);
% np = normalizenpmat4(reorder_np_cell_to_array(np));
% 
% [errR_mean, errR_median, errT_mean, errT_median] = comparison_for_calibrated_cameras(reorder_np_array_to_cell(cameras), reorder_np_array_to_cell(np));
% res.errR_mean = errR_mean;
% res.errR_median = errR_median;
% res.errT_mean = errT_mean / t_scale;
% res.errT_median = errT_median / t_scale;
% res.res2 = res2;
% 
% cameras_est = reorder_np_array_to_cell(cameras);
% end
% 
% function npmat = normalizenpmat4(npmat)
%     n = size(npmat,1)/3;
% %     for i = 1:n
% %        npmat(3*(i-1)+1:3*i,4) = npmat(3*(i-1)+1:3*i,4) / norm(npmat(3*(i-1)+1:3*i,4),'fro'); 
% %     end
%     for i = 1:n
%        npmat(3*(i-1)+1:3*i,:) = npmat(3*(i-1)+1:3*i,:) / norm(npmat(3*(i-1)+1:3*i,1:3),'fro'); 
%     end
% end
% 
% function [metrics, metrics_rotation, metrics_translation, good_measurements] = compare_block_tensor(cTest, n, cameras_gt,t_scale)
% 
% % T_gt = generate_block_trifocal_tensor(cameras_gt, n);
% 
% metrics = zeros(n,n,n);
% metrics_rotation = zeros(n,n,n);
% metrics_translation = zeros(n,n,n);
% good_measurements = zeros(n,n,n);
% parfor i= 1:n
%     for j = 1:n
%         for k = 1:n
%             if nnz(cTest{i,j,k}) > 0
%                 T_gt = T_from_P({cameras_gt{k}, cameras_gt{i}, cameras_gt{j}});
% 
%                 cur_rel_err = scaled_diff_frobenius_norm(cTest{i,j,k}, T_gt);
%                 %%% Try filtering the dataset
%                 % if cur_rel_err > 0.2
%                 %     continue;
%                 % end
% 
%                 % good_measurements(i,j,k) = 1;
%                 metrics(i,j,k) = cur_rel_err;
%                 [errR_mean, errT_mean] = angle_diff(cTest{i,j,k}, T_gt);
%                 metrics_rotation(i,j,k) = errR_mean;
%                 metrics_translation(i,j,k) = errT_mean / t_scale;
%                 
%             end
%         end
%     end
%     fprintf("iteration %d\n", i);
% end
% 
% 
% end