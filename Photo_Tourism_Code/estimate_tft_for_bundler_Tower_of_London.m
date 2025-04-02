% estimate_tft_for_bundler
% addpath(genpath(pwd))
dataset = "Tower_of_London";
load("External Packages and Toolboxes/SfM_CVPR2017_code/Data/"+dataset+"_data.mat")
Jmat = diag([1 1 -1]);
n = data.n;
Rmat_orig = zeros(3,3,n);
tmat_orig = zeros(3,1,n);
ccams = zeros(3,4,n);
cameras_gt = cell(1,n);
t_scale = sqrt(n) / frob(data.t);
for i = 1:n
    Rmat_orig(:,:,i) = data.R(:,:,i)';
    tmat_orig(:,:,i) = - Rmat_orig(:,:,i) * data.t(:,i) * t_scale;
    ccams(:,:,i) = Rmat_orig(:,:,i)' * [eye(3), -tmat_orig(:,:,i)];
    cameras_gt{i} = ccams(:,:,i);
end

good_index = [];
for i = 1:n
    Rmat_orig(:,:,i) = data.R(:,:,i)';
    tmat_orig(:,:,i) = - Rmat_orig(:,:,i) * data.t(:,i) * t_scale;
    ccams(:,:,i) = Rmat_orig(:,:,i)' * [eye(3), -tmat_orig(:,:,i)];
    cameras_gt{i} = ccams(:,:,i);
    if nnz(cameras_gt{i}) > 0
        good_index = [good_index, i];
    end
end

save("Photo_Tourism_cleaning/" + dataset+"/"+dataset+"good_index.mat", "good_index")

%% filter the void estimates
E_est = data.E_est(good_index, good_index);
cameras_gt = cameras_gt(good_index);
ccams = ccams(:,:,good_index);
Rmat_orig = Rmat_orig(:,:,good_index);
tmat_orig = tmat_orig(:,:,good_index);
n = length(cameras_gt);

fprintf("Just cleared all the void initial estimates\n")
clear data;

%% calculate cameras from essential matrices
tft_estimated_cams = cell(n,n,n);
data_count = zeros(n,n,n);
parfor i= 1:n
    for j = 1:n
        for k = 1:n
            if i < j && j < k && size(E_est{j,i},1) > 0 && size(E_est{k,j},1) > 0 && size(E_est{k,i},1) > 0

                [P1,P2,P3] = tft_from_3_essmat(E_est{j,i}, E_est{k,i}, E_est{k,j});
                tft_estimated_cams{i,j,k} = {P1,P2,P3};
                data_count(i,j,k) = 1;
            end
        end
    end
    fprintf("iteration %d \n", i);
end

for i = 1:n
    for j = 1:n
        for k = 1:n
            if i<j && j < k && data_count(i,j,k)
                data_count(j,j,k) = 1;
                data_count(k,k,j) = 1;
                data_count(i,i,j) = 1;
                data_count(j,j,i) = 1;
                data_count(i,i,k) = 1;
                data_count(k,k,i) = 1;

                data_count(j,k,i) = 1;
                data_count(k,j,i) = 1;
                data_count(i,k,j) = 1;
                data_count(j,i,j) = 1;
                data_count(j,i,k) = 1;
            end
        end
    end
end

%% We first truncate the size of cTest then calculate it 
datacounts = sum(sum(data_count));
first_filtering = datacounts > (0.03 * n^2);
first_filtering = find(first_filtering);
if length(first_filtering) > 225
    [~,first_filtering]  = maxk(datacounts, 225);
end
save("Photo_Tourism_cleaning/" + dataset+"/"+dataset+"first_filtering_indices.mat", "first_filtering")

fprintf("The size of the first filtering is %d\n", length(first_filtering));
tft_estimated_cams = tft_estimated_cams(first_filtering, first_filtering, first_filtering);
cameras_gt = cameras_gt(first_filtering);
E_est = E_est(first_filtering, first_filtering);
n = length(first_filtering);
data_count = data_count(first_filtering, first_filtering, first_filtering);

cTest = cell(n,n,n);
for i = 1:n
    for j = 1:n
        for k = 1:n
            if i <j && j < k && size(tft_estimated_cams{i,j,k}, 1)>0
                P1 = tft_estimated_cams{i,j,k}{1};
                P2 = tft_estimated_cams{i,j,k}{2};
                P3 = tft_estimated_cams{i,j,k}{3};


                cTest{j,k,i} = T_from_P({P1,P2,P3});
                cTest{k,j,i} = T_from_P({P1,P3,P2});
                cTest{i,k,j} = T_from_P({P2,P1,P3});
                cTest{k,i,j} = T_from_P({P2,P3,P1});
                cTest{i,j,k} = T_from_P({P3,P1,P2});
                cTest{j,i,k} = T_from_P({P3,P2,P1});

            end
        end
    end
end
clear tft_estimated_cams datacounts data_count;
% save("Photo_Tourism_cleaning/"+dataset+"/"+dataset+"_tft_from_essmat.mat", "cTest", '-v7.3');

fprintf("loading the estimated trifocal tensors \n");
% load("Photo_Tourism_cleaning/"+dataset+"/"+dataset+"_tft_from_essmat.mat")
parfor i = 1:n
    for j = 1:n
        if nnz(E_est{i,j}) > 0
            cTest{j,j,i} = reorder_iij_from_fundamental(E_est{i,j}); 
        end
    end
end


diary Photo_Tourism_cleaning/Tower_of_London/output.log
fprintf("Total number of cameras is %d \n", length(first_filtering));
try
    fprintf("%s starting to compare trifocal tensors \n", dataset);
catch
end


T = final_T_cell_to_array(cTest);
clear cTest;
fprintf("Starting to synchronize\n");
tic
final_T = heuristic_hosvd_imputation(T, cameras_gt, t_scale);
toc
save("Photo_Tourism_cleaning/"+dataset+"/"+dataset+"_final_tensor_filtered.mat","cameras_gt", "final_T", '-v7.3')


% load("EPFL_code/"+dataset+"_final_tensor.mat")

[res, cameras_est] = residual_hosvd(final_T, cameras_gt, t_scale);
disp('HOSVD');
disp('Relative Cam Difference | Mean Location Error | Median Location Error | Mean Rotation Error | Median Rotation Error |');
fprintf('%.4f %.4f %.4f %.4f %.4f \n', res.res2, res.errT_mean, res.errT_median, res.errR_mean, res.errR_median);
diary off

save("Photo_Tourism_cleaning/"+dataset+"/"+dataset+"_cameras_est_filtered.mat", "res", "cameras_est", '-v7.3');


function cams_good = get_cams_from_final_T(final_T, good_index)
[U,S,sv] = mlsvd(final_T);
cs = U{1}(:,1:4);
B = U{2}(:,1:4);
C = U{3}(:,1:6);
G = S(1:4, 1:4, 1:6);
cur_t2 = tmprod(G, {cs,B,C}, 1:3);


cur_t2(isnan(cur_t2)) = 0;
cameras = retrieve_cameras_hosvd(cur_t2);

cams_good = zeros(3*length(good_index), 4);
for i = 1:length(good_index)
   cams_good(3*(i-1)+1:3*i,:) = cameras(3*(i-1)+1:3*i, :); 
end

end

function [res, cameras_est] = residual_hosvd(cur_t, np, t_scale)
[U,S,sv] = mlsvd(cur_t);
cs = U{1}(:,1:4);
B = U{2}(:,1:4);
C = U{3}(:,1:6);
G = S(1:4, 1:4, 1:6);
cur_t2 = tmprod(G, {cs,B,C}, 1:3);


cur_t2(isnan(cur_t2)) = 0;
cameras = retrieve_cameras_hosvd(cur_t2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uncali = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[dd2,nn2,Hf, npmat1, npmat2] = compare_projection_matrices(reorder_np_array_to_cell(cameras), np, uncali);
res2 = dd2 / nn2;

cameras = normalizenpmat4(cameras);
np = normalizenpmat4(reorder_np_cell_to_array(np));

[errR_mean, errR_median, errT_mean, errT_median] = comparison_for_calibrated_cameras(reorder_np_array_to_cell(cameras), reorder_np_array_to_cell(np));
res.errR_mean = errR_mean;
res.errR_median = errR_median;
res.errT_mean = errT_mean / t_scale;
res.errT_median = errT_median / t_scale;
res.res2 = res2;

cameras_est = reorder_np_array_to_cell(cameras);
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

function [metrics, metrics_rotation, metrics_translation, good_measurements] = compare_block_tensor(cTest, n, cameras_gt,t_scale)

% T_gt = generate_block_trifocal_tensor(cameras_gt, n);

metrics = zeros(n,n,n);
metrics_rotation = zeros(n,n,n);
metrics_translation = zeros(n,n,n);
good_measurements = zeros(n,n,n);
parfor i= 1:n
    for j = 1:n
        for k = 1:n
            if nnz(cTest{i,j,k}) > 0
                T_gt = T_from_P({cameras_gt{k}, cameras_gt{i}, cameras_gt{j}});

                cur_rel_err = scaled_diff_frobenius_norm(cTest{i,j,k}, T_gt);
                %%% Try filtering the dataset
                % if cur_rel_err > 0.2
                %     continue;
                % end

                % good_measurements(i,j,k) = 1;
                metrics(i,j,k) = cur_rel_err;
                [errR_mean, errT_mean] = angle_diff(cTest{i,j,k}, T_gt);
                metrics_rotation(i,j,k) = errR_mean;
                metrics_translation(i,j,k) = errT_mean / t_scale;
                
            end
        end
    end
    fprintf("iteration %d\n", i);
end


end



% % % estimate_tft_for_bundler
% % % addpath(genpath(pwd))
% % dataset = "Yorkminster";
% % load("External Packages and Toolboxes/SfM_CVPR2017_code/Data/"+dataset+"_data.mat")
% % Jmat = diag([1 1 -1]);
% % n = data.n;
% % Rmat_orig = zeros(3,3,n);
% % tmat_orig = zeros(3,1,n);
% % ccams = zeros(3,4,n);
% % cameras_gt = cell(1,n);
% % t_scale = sqrt(n) / frob(data.t);
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
% % 
% % cameras_gt = cameras_gt(good_index);
% % ccams = ccams(:,:,good_index);
% % Rmat_orig = Rmat_orig(:,:,good_index);
% % tmat_orig = tmat_orig(:,:,good_index);
% % n = length(cameras_gt);
% % save("Photo_Tourism_cleaning/Yorkminster/Yorkminster_initial_good_index.mat", "good_index", "cameras_gt");
% % 
% % 
% % % tft_estimated_cams = cell(n,n,n);
% % % parfor i= 1:n
% % %     for j = 1:n
% % %         for k = 1:n
% % %             if i < j && j < k && size(data.E_est{j,i},1) > 0 && size(data.E_est{k,j},1) > 0 && size(data.E_est{k,i},1) > 0
% % % 
% % %                 [P1,P2,P3] = tft_from_3_essmat(data.E_est{j,i}, data.E_est{k,i}, data.E_est{k,j});
% % %                 tft_estimated_cams{i,j,k} = {P1,P2,P3};
% % % 
% % %             end
% % %         end
% % %     end
% % %     fprintf("iteration %d \n", i);
% % % end
% % % cTest = cell(n,n,n);
% % % for i = 1:n
% % %     for j = 1:n
% % %         for k = 1:n
% % %             if i <j && j < k && size(tft_estimated_cams{i,j,k}, 1)>0
% % %                 P1 = tft_estimated_cams{i,j,k}{1};
% % %                 P2 = tft_estimated_cams{i,j,k}{2};
% % %                 P3 = tft_estimated_cams{i,j,k}{3};
% % % 
% % % 
% % %                 cTest{j,k,i} = T_from_P({P1,P2,P3});
% % %                 cTest{k,j,i} = T_from_P({P1,P3,P2});
% % %                 cTest{i,k,j} = T_from_P({P2,P1,P3});
% % %                 cTest{k,i,j} = T_from_P({P2,P3,P1});
% % %                 cTest{i,j,k} = T_from_P({P3,P1,P2});
% % %                 cTest{j,i,k} = T_from_P({P3,P2,P1});
% % % 
% % %             end
% % %         end
% % %     end
% % % end
% % % clear tft_estimated_cams;
% % % save("Photo_Tourism_cleaning/"+dataset+"/"+dataset+"_tft_from_essmat.mat", "cTest", '-v7.3');
% % 
% % fprintf("loading the estimated trifocal tensors \n");
% % load("Photo_Tourism_cleaning/"+dataset+"/"+dataset+"_tft_from_essmat.mat");
% % 
% % 
% % parfor i = 1:n
% %     for j = 1:n
% %         if nnz(data.E_est{i,j}) > 0
% %             cTest{j,j,i} = reorder_iij_from_fundamental(data.E_est{i,j}); 
% %         end
% %     end
% % end
% % 
% % %% We first truncate the size of cTest
% % first_filtering = zeros(1,n);
% % parfor i = 1:n
% %     count = 0;
% %     for j = 1:n
% %         for k = 1:n
% %             if size(cTest{i,j,k},1) > 0
% %                 count = count + 1;
% %             end
% % 
% %         end
% %     end
% %     if count >= 0.1 * n^2
% %         first_filtering(i) = count;
% %     end
% % end
% % 
% % first_filtering = find(first_filtering);
% % if length(first_filtering) > 300
% %     [~,first_filtering]  = maxk(first_filtering, 300);
% % end
% % fprintf("The size of the filtering is %d \n", length(first_filtering));
% % cTest = cTest(first_filtering, first_filtering, first_filtering);
% % cameras_gt = cameras_gt(first_filtering);
% % 
% % diary Photo_Tourism_cleaning/Yorkminster/output.log
% % try
% %     fprintf("%s starting to compare trifocal tensors \n", dataset);
% % catch
% % end
% % % [metrics, metrics_rotation, metrics_translation, good_measurements] = compare_block_tensor(cTest, n, cameras_gt, t_scale);
% % % figure;
% % % subplot(1,3,1)
% % % histogram(metrics(metrics >0),50)
% % % title("Method Relative TFT")
% % % 
% % % subplot(1,3,2)
% % % metrics_rotation = metrics_rotation(metrics_rotation > 0);
% % % histogram(metrics_rotation(:),50)
% % % title("Method Rotation Errors")
% % % 
% % % subplot(1,3,3)
% % % metrics_translation = metrics_translation(metrics_translation > 0);
% % % histogram(metrics_translation(:),50)
% % % title("Method Translation Errors")
% % 
% % T = final_T_cell_to_array(cTest);
% % clear cTest;
% % fprintf("Starting to synchronize\n");
% tic
% final_T = heuristic_hosvd_imputation(T, cameras_gt, t_scale);
% toc
% save("Photo_Tourism_cleaning/"+dataset+"/"+dataset+"_final_tensor_filtered.mat","cameras_gt", "final_T", "good_index", '-v7.3')
% 
% 
% [res, cameras_est] = residual_hosvd(final_T, cameras_gt, t_scale);
% disp('HOSVD');
% disp('Relative Cam Difference | Mean Location Error | Median Location Error | Mean Rotation Error | Median Rotation Error |');
% fprintf('%.4f %.4f %.4f %.4f %.4f \n', res.res2, res.errT_mean, res.errT_median, res.errR_mean, res.errR_median);
% diary off
% 
% save("Photo_Tourism_cleaning/"+dataset+"/"+dataset+"_cameras_est_filtered.mat", 'res', 'cameras_est', '-v7.3');
% 
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
% end
