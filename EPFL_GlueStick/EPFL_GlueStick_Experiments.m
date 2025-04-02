path_to_data=strcat('EPFL_code/Data/',dataset,'/',dataset,'_dense/');

initial_sample_size=100;
bundle_adj_size=50;
repr_err_th=1;



% load("EPFL_code/EPFL_"+dataset+"_FCC_matched_full.mat")
% dataset = "FountainP11";
% dataset_name = "Herz_P8";
% dataset_name = "Entry_P10";
% dataset_name = "Castle_P19";
dataset_name = "Fountain_P11";
n = 11;

load("EPFL_GlueStick/"+dataset+"_two_view_data.mat");

[triplet_points, triplet_lines, n] = EPFL_process_data(dataset_name, n);

corresp_by_triplet = triplet_points;
% corresp_by_triplet = matched_triplets_epfl;

n = size(corresp_by_triplet,1);
for i = 1:n
    for j = 1:n
        for k = 1:n
            if size(corresp_by_triplet{i,j,k},1) > 0
                corresp_by_triplet{i,j,k} = unique(corresp_by_triplet{i,j,k}, 'rows');
            end
        end
    end
end

im_names

[Ks, Rs, ts, ccams, cTgt] = EPFL_ground_truth(im_names, path_to_data, t_scale);
bad_reprojection = [];
%% methods to evaluate

methods={@STETFTPoseEstimation};     

methods_to_test = [1];

%% error vectors
repr_err = zeros(n^3,length(methods),2);
rot_err  = zeros(n^3,length(methods),2);
t_err    = zeros(n^3,length(methods),2);
iter     = zeros(n^3,length(methods),2);
time     = zeros(n^3,length(methods),2);
time_estimate = zeros(n^3,4);

%% evaluation
cTest = cell(n,n,n);
it = 1;
tic
for i = 1:n
    for j = 1:n
        for k = 1:n
            if i < j && j < k && (size(data.E_est{i,j},1) > 0) && ...
                    (size(data.E_est{i,k},1) > 0) && (size(data.E_est{j,k},1) > 0)

            % if i < j && j < k && (size(corresp_by_triplet{i,j,k},1) > 10)

                % Triplet information and correspondances

                im1 = i; im2 = j; im3 = k; 
                Corresp=corresp_by_triplet{im1,im2,im3}';

                if size(Corresp,2) < 10
                    continue
                end

                N=size(Corresp,2);
                fprintf('Triplet %d/%d (%d,%d,%d) with %d matching points.\n',...
                    it,n^3,im1,im2,im3,N);

                % Ground truth poses and calibration
                [K1,R1_true,t1_true,im_size]=readCalibrationOrientation_EPFL(path_to_data,im_names{im1});
                [K2,R2_true,t2_true]=readCalibrationOrientation_EPFL(path_to_data,im_names{im2});
                [K3,R3_true,t3_true]=readCalibrationOrientation_EPFL(path_to_data,im_names{im3});
                CalM=[K1;K2;K3];
                R_t0={[R2_true*R1_true.', t2_true-R2_true*R1_true.'*t1_true],...
                      [R3_true*R1_true.', t3_true-R3_true*R1_true.'*t1_true]};
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Corresp_inliers = Corresp;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if size(Corresp,2) > 1500
                    rng(it);
                    init_sample=randsample(1:N,1500);
                    rng(it);
                    Corresp_init=Corresp_inliers(:,init_sample);
                    Corresp = Corresp_init;
                end

                % estimate the pose using all methods
                fprintf('method... ');
                for m=methods_to_test
                    fprintf('%d ',m);

                    % if there are not enough matches for initial estimation, inf error
                    if (m>6 && N<8) || N<7 
                        repr_err(it,m,:)=inf;    rot_err(it,m,:)=inf;
                        t_err(it,m,:)=inf;       iter(it,m,:)=inf;
                        time(it,m,:)=inf;
                        continue;
                    end

                    % % STE POSE ESTIMATION      TIMING 1
                    t0=cputime;
                    [R_t_2,R_t_3,Reconst,T,nit, distances]=methods{m}(Corresp,CalM);   % STE to find 
                    t=cputime-t0;
                    time_estimate(it, 1) = t;

                    t0 = cputime;
                    [T, distances2, Corresp_inliers, R_t_2, R_t_3] = STE_local_optimization(T, Corresp', distances, 0.4, CalM);
                    t=cputime-t0;
                    time_estimate(it, 2) = t;
                    % % STE POSE ESTIMATION       

                    Corresp_inliers = Corresp_inliers';
                    N = size(Corresp_inliers,2);

                    rng(it);
                    init_sample=randsample(1:N,min(initial_sample_size,N));
                    rng(it);
                    ref_sample=randsample(init_sample,min(bundle_adj_size,length(init_sample)));
                    Corresp_inliers=Corresp_inliers(:,init_sample);
                    Corresp_ref=Corresp_inliers(:,ref_sample);


                    % reprojection error with all inliers
%                     repr_err(it,m,1)= ReprError({CalM(1:3,:)*eye(3,4),...
%                         CalM(4:6,:)*R_t_2,CalM(7:9,:)*R_t_3},Corresp_inliers);
%                     % angular errors
%                     [rot2_err,t2_err]=AngError(R_t0{1},R_t_2);
%                     [rot3_err,t3_err]=AngError(R_t0{2},R_t_3);
%                     rot_err(it,m,1)=(rot2_err+rot3_err)/2;
%                     t_err(it,m,1)=(t2_err+t3_err)/2;
%                     % iterations & time
%                     iter(it,m,1)=nit; time(it,m,1)=t;


                    % % Apply Bundle Adjustment TIMING 2
                    fprintf('(ref)... ');
                    t0=cputime;
                    [R_t_ref,~,nit,repr_errBA]=BundleAdjustment(CalM,...
                        [eye(3,4);R_t_2;R_t_3],Corresp_ref);
                    t=cputime-t0;
                    time_estimate(it, 3) = t;

                    % reprojection error with all inliers TIMING 3
                    t0 = cputime;
                    repr_err(it,m,2)= ReprError({CalM(1:3,:)*R_t_ref(1:3,:),...
                        CalM(4:6,:)*R_t_ref(4:6,:),...
                        CalM(7:9,:)*R_t_ref(7:9,:)},Corresp_inliers);

                    ref_R_t = {R_t_ref(1:3,:),R_t_ref(4:6,:), R_t_ref(7:9,:)};
                    if repr_err(it,m,2) < 1
                        cTest{im2,im3,im1} = T_from_P({R_t_ref(1:3,:),R_t_ref(4:6,:), R_t_ref(7:9,:)});
                        cTest{im3,im2,im1} = T_from_P({R_t_ref(1:3,:),R_t_ref(7:9,:), R_t_ref(4:6,:)});
                        cTest{im1,im3,im2} = T_from_P({R_t_ref(4:6,:),R_t_ref(1:3,:), R_t_ref(7:9,:)});
                        cTest{im3,im1,im2} = T_from_P({R_t_ref(4:6,:),R_t_ref(7:9,:), R_t_ref(1:3,:)});
                        cTest{im1,im2,im3} = T_from_P({R_t_ref(7:9,:),R_t_ref(1:3,:), R_t_ref(4:6,:)});
                        cTest{im2,im1,im3} = T_from_P({R_t_ref(7:9,:),R_t_ref(4:6,:), R_t_ref(1:3,:)});
                    else
                        bad_reprojection = [bad_reprojection; i,j,k];
                    end
                    t = cputime - t0;
                    time_estimate(it, 4) = t;


                    % angular errors
%                     [rot2_err,t2_err]=AngError(R_t0{1},R_t_ref(4:6,:));
%                     [rot3_err,t3_err]=AngError(R_t0{2},R_t_ref(7:9,:));
%                     rot_err(it,m,2)=(rot2_err+rot3_err)/2;
%                     t_err(it,m,2)=(t2_err+t3_err)/2;
%                     % iterations & time
%                     iter(it,m,2)=nit; time(it,m,2)=t;

                end
                it = it + 1; 
            end
        end
    end
end
toc

% load("EPFL_GlueStick/"+dataset+"_two_view_data.mat");
% for i = 1:n
%     for j = 1:n
%         if i ~= j
%             cTest{i,i,j} = reorder_iij_from_fundamental(data.E_est{i,j});
%         end
%     end
% end


figure;
subplot(2,2,1)
histogram(time_estimate(time_estimate(:,1) >0,1),50)
title("STE estimate inlier time")

subplot(2,2,2)
histogram(time_estimate(time_estimate(:,2) >0,2),50)
title("linear estimator time")

subplot(2,2,3)
histogram(time_estimate(time_estimate(:,3) >0,3),50)
title("Bundle Adjustment Time")

subplot(2,2,4)
histogram(time_estimate(time_estimate(:,4) >0,4),50)
title("filter and calculate time")




cameras_gt = cell(1,n);
for i = 1:n
   cameras_gt{i} = ccams(:,:,i); 
end

T = final_T_cell_to_array(cTest);
tic
final_T = heuristic_hosvd_imputation(T, cameras_gt, t_scale);
toc
[res, SO3Mats_corr2, camorigR, t_fit2, torig_centered] = residual_hosvd(final_T, cameras_gt, t_scale);
disp('HO-rSTE');
disp('Relative Cam Difference | Mean Location Error | Median Location Error | Mean Rotation Error | Median Rotation Error |');
fprintf('%.4f %.4f %.4f %.4f %.4f \n', res.res2, res.errT_mean, res.errT_median, res.errR_mean, res.errR_median);
save(dataset+"_trifocal_tensor_data.mat")

% all_filter_indices = [1:5;4:8;7:11;10:14;13:17;16:19,1];
% 
% calced_cams = {};
% 
% for i = 1:size(all_filter_indices,1)
%     filter_indices = all_filter_indices(i,:);
%     cur_Test = cTest(filter_indices, filter_indices, filter_indices);
%     cur_cams_gt = cameras_gt(filter_indices);
% 
%     tic
%     final_T = heuristic_hosvd_imputation(final_T_cell_to_array(cur_Test), cur_cams_gt, t_scale);
%     toc
%     % save("EPFL_code/"+dataset+"_final_tensor.mat", "T", "cameras_gt", "final_T")
%     calced_cams{i} = get_calced_cams(final_T);
%     % load("EPFL_code/"+dataset+"_final_tensor.mat")
%     [res, SO3Mats_corr2, camorigR, t_fit2, torig_centered] = residual_hosvd(final_T, cur_cams_gt, t_scale);
%     disp('HO-rSTE');
%     disp('Relative Cam Difference | Mean Location Error | Median Location Error | Mean Rotation Error | Median Rotation Error |');
%     fprintf('%.4f %.4f %.4f %.4f %.4f \n', res.res2, res.errT_mean, res.errT_median, res.errR_mean, res.errR_median);
% end
% 
% final_cams = synchronize_all_subproblems(calced_cams);

% cameras = normalizenpmat4(final_cams);
% np = normalizenpmat4(reorder_np_cell_to_array(cameras_gt));
% 
% cameras = cameras(1:57,:);
% 
% [errR_mean, errR_median, errT_mean, errT_median, SO3Mats_corr2, camorigR, t_fit2, torig_centered] = comparison_for_calibrated_cameras(reorder_np_array_to_cell(cameras), reorder_np_array_to_cell(np));
% res.errR_mean = errR_mean;
% res.errR_median = errR_median;
% res.errT_mean = errT_mean / t_scale;
% res.errT_median = errT_median / t_scale;
% disp('HO-rSTE');
% disp('Relative Cam Difference | Mean Location Error | Median Location Error | Mean Rotation Error | Median Rotation Error |');
% fprintf('%.4f %.4f %.4f %.4f %.4f \n', 0, res.errT_mean, res.errT_median, res.errR_mean, res.errR_median);

% save(dataset+"_cameras_est.mat")
% save(dataset+"_trifocal_tensor_data.mat")

function final_cams = synchronize_all_subproblems(calced_cams)

final_cams = calced_cams{1};
for i = 1:length(calced_cams)
    if i < length(calced_cams)
        [Hf ,Psa, vscales] = findHomography_large({calced_cams{i+1}(1:3,:), calced_cams{i+1}(4:6,:)}, ...
            {calced_cams{i}(10:12,:), calced_cams{i}(13:15,:)});
        calced_cams{i+1} = calced_cams{i+1} * Hf;
        final_cams = [final_cams; calced_cams{i+1}(7:end,:)];
    end
end




end



function cameras = get_calced_cams(cur_t)
[U,S,sv] = mlsvd(cur_t);
cs = U{1}(:,1:4);
B = U{2}(:,1:4);
C = U{3}(:,1:6);
G = S(1:4, 1:4, 1:6);
cur_t2 = tmprod(G, {cs,B,C}, 1:3);


cur_t2(isnan(cur_t2)) = 0;
cameras = retrieve_cameras_hosvd(cur_t2);

end


function [res, SO3Mats_corr2, camorigR, t_fit2, torig_centered] = residual_hosvd(cur_t, np, t_scale)
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

[errR_mean, errR_median, errT_mean, errT_median, SO3Mats_corr2, camorigR, t_fit2, torig_centered] = comparison_for_calibrated_cameras(reorder_np_array_to_cell(cameras), reorder_np_array_to_cell(np));
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


