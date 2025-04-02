% % % Script adopted from submission "A Critical Review of the Trifocal Tensor Estimation"
% % 
% % % Copyright (c) 2017 Laura F. Julia <laura.fernandez-julia@enpc.fr>
% % % All rights reserved.
% % %
% % % This program is free software: you can redistribute it and/or modify
% % % it under the terms of the GNU General Public License as published by
% % % the Free Software Foundation, either version 3 of the License, or
% % % (at your option) any later version.
% % % 
% % % This program is distributed in the hope that it will be useful,
% % % but WITHOUT ANY WARRANTY; without even the implied warranty of
% % % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% % % GNU General Public License for more details.
% % % 
% % % You should have received a copy of the GNU General Public License
% % % along with this program.  If not, see <http://www.gnu.org/licenses/>.

% % %% Some parameters
 
path_to_data=strcat('EPFL_code/Data/',dataset,'/',dataset,'_dense/');

initial_sample_size=100;
bundle_adj_size=50;
repr_err_th=1;


% Recover correspondances

% % % % % corresp_file = matfile(strcat(path_to_data,'Corresp_triplets','.mat'));
% % % % % indexes_sorted = corresp_file.indexes_sorted;
% % % % % corresp_by_triplet = corresp_file.Corresp;
% % % % % im_names = corresp_file.im_names;
% % % % % clear corresp_file;

load("EPFL_code/EPFL_"+dataset+"_FCC_matched_full.mat")
corresp_by_triplet = matched_triplets_epfl;
 
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

[Ks, Rs, ts, ccams, cTgt] = EPFL_ground_truth(im_names, path_to_data, t_scale);
bad_reprojection = [];
%% methods to evaluate

methods={...
%     @LinearTFTPoseEstimation,...    % 1 - TFT - Linear estimation
        @STETFTPoseEstimation
%     @ResslTFTPoseEstimation     % 2 - TFT - Ressl
%     @NordbergTFTPoseEstimation,...  % 3 - TFT - Nordberg
%     @FaugPapaTFTPoseEstimation,...  % 4 - TFT - Faugeras&Papadopoulo 
%     @PiPoseEstimation,...           % 5 - Pi matrices - Ponce&Hebert
%     @PiColPoseEstimation,...        % 6 - Pi matrices - Ponce&Hebert for collinear cameras
%     @LinearFPoseEstimation,...      % 7 - Fundamental matrices - Linear estimation
%     @OptimFPoseEstimation
    };         % 8 - Fundamental matrices - Optimized

% methods_to_test=[1:5,7:8]; % no method for collinear cameras
methods_to_test = [1];

%% error vectors
repr_err = zeros(n^3,length(methods),2);
rot_err  = zeros(n^3,length(methods),2);
t_err    = zeros(n^3,length(methods),2);
iter     = zeros(n^3,length(methods),2);
time     = zeros(n^3,length(methods),2);

%% evaluation
cTest = cell(n,n,n);
it = 1;
for i = 1:n
    for j = 1:n
        for k = 1:n
            if i < j && j < k && (size(corresp_by_triplet{i,j,k},1) > 10)

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

            % %     % Discart correspondances with repr_err > 1 pix
            %     Reconst0=triangulation3D({K1*eye(3,4),K2*R_t0{1},K3*R_t0{2}},Corresp);
            %     Reconst0=bsxfun(@rdivide,Reconst0(1:3,:),Reconst0(4,:));
            %     Corresp_new=project3Dpoints(Reconst0,{K1*eye(3,4),K2*R_t0{1},K3*R_t0{2}});
            %     residuals=Corresp_new-Corresp; % reprojection error
            %     Corresp_inliers=Corresp(:, sum(abs(residuals)>repr_err_th,1)==0 );
            %     fprintf("%f percentage of inliers\n", size(Corresp_inliers,2) /N);
            %     N=size(Corresp_inliers,2);
            %     REr=ReprError({K1*eye(3,4),K2*R_t0{1},K3*R_t0{2}},Corresp_inliers);
            %     fprintf('%d valid correspondances with reprojection error %f.\n',N,REr);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Corresp_inliers = Corresp;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if size(Corresp,2) > 1500
                    rng(it);
                    init_sample=randsample(1:N,1500);
                    rng(it);
            %         ref_sample=randsample(init_sample,min(bundle_adj_size,length(init_sample)));
                    Corresp_init=Corresp_inliers(:,init_sample);
            %         Corresp_ref=Corresp_inliers(:,ref_sample);
                    Corresp = Corresp_init;
                end

                % Samples for initial estimation and bundle adjustment
            %     rng(it);
            %     init_sample=randsample(1:N,min(initial_sample_size,N));
            %     rng(it);
            %     ref_sample=randsample(init_sample,min(bundle_adj_size,length(init_sample)));
            %     Corresp_init=Corresp_inliers(:,init_sample);
            %     Corresp_ref=Corresp_inliers(:,ref_sample);
            %    init_sample=randsample(1:N,min(initial_sample_size,N));

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

                    % % Pose estimation by method m, measuring time            
                    t0=cputime;
                    [R_t_2,R_t_3,Reconst,T,nit, distances]=methods{m}(Corresp,CalM);   % STE to find 
                    [T, distances2, Corresp_inliers, R_t_2, R_t_3] = STE_local_optimization(T, Corresp', distances, 0.4, CalM);

                    t=cputime-t0;

                    Corresp_inliers = Corresp_inliers';
                    N = size(Corresp_inliers,2);

                    rng(it);
                    init_sample=randsample(1:N,min(initial_sample_size,N));
                    rng(it);
                    ref_sample=randsample(init_sample,min(bundle_adj_size,length(init_sample)));
                    Corresp_inliers=Corresp_inliers(:,init_sample);
                    Corresp_ref=Corresp_inliers(:,ref_sample);


                    % reprojection error with all inliers
                    repr_err(it,m,1)= ReprError({CalM(1:3,:)*eye(3,4),...
                        CalM(4:6,:)*R_t_2,CalM(7:9,:)*R_t_3},Corresp_inliers);
                    % angular errors
                    [rot2_err,t2_err]=AngError(R_t0{1},R_t_2);
                    [rot3_err,t3_err]=AngError(R_t0{2},R_t_3);
                    rot_err(it,m,1)=(rot2_err+rot3_err)/2;
                    t_err(it,m,1)=(t2_err+t3_err)/2;
                    % iterations & time
                    iter(it,m,1)=nit; time(it,m,1)=t;


                    % % Apply Bundle Adjustment
                    fprintf('(ref)... ');
                    t0=cputime;
                    [R_t_ref,~,nit,repr_errBA]=BundleAdjustment(CalM,...
                        [eye(3,4);R_t_2;R_t_3],Corresp_ref);
                    t=cputime-t0;

                    % reprojection error with all inliers
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

                    % angular errors
                    [rot2_err,t2_err]=AngError(R_t0{1},R_t_ref(4:6,:));
                    [rot3_err,t3_err]=AngError(R_t0{2},R_t_ref(7:9,:));
                    rot_err(it,m,2)=(rot2_err+rot3_err)/2;
                    t_err(it,m,2)=(t2_err+t3_err)/2;
                    % iterations & time
                    iter(it,m,2)=nit; time(it,m,2)=t;

                end
                it = it + 1; 
            end
        end
    end
end

% % load("EPFL_code/EPFL_fountainp11_debug_comparison.mat");
% save("EPFL_code/"+dataset+"_cTest_triples.mat","cTest");
% 
% load("EPFL_code/"+dataset+"_cTest_triples.mat");
% 

cameras_gt = cell(1,n);
for i = 1:n
   cameras_gt{i} = ccams(:,:,i); 
end
% 
% [metrics, metrics_rotation, metrics_translation] = compare_block_tensor(cTest, n, cameras_gt);
% figure;
% subplot(1,3,1)
% histogram(metrics(metrics >0),50)
% title("Method Relative TFT")
% 
% subplot(1,3,2)
% metrics_rotation = metrics_rotation(metrics_rotation > 0);
% histogram(metrics_rotation(:),50)
% title("Method Rotation Errors")
% 
% subplot(1,3,3)
% metrics_translation = metrics_translation(metrics_translation > 0);
% histogram(metrics_translation(:),50)
% title("Method Translation Errors")



T = final_T_cell_to_array(cTest);

%% HARD CODE FOR HERZP25, AS IT IS THE ONLY ONE THAT HAS CAMERAS WITH
%% MISSING MEASUREMENTS
if strcmp(dataset, "HerzP25")
    good_indices = [];
    for i = 1:(3*n)
        if nnz(T(i,:,:)) == 0 && ~mod(i,3)
            continue;
        else
            if nnz(T(i,:,:))>0  && ~mod(i,3)
                good_indices = [good_indices,i/3];
            end
        end
    end
%     good_indices = [1,2,3,4,5,6,7,8,9,10,11,12,13,15,16,17,18,19,20,21,22,23,24,25];
    cTest2 = cTest(good_indices, good_indices, good_indices);
    save("EPFL_code/HerzP25_good_index.mat","good_indices");
    cameras_gt2 = cameras_gt(good_indices);
    fprintf("We actually have only %d cameras out of 25 with estimates\n", length(good_indices));
    T2 = final_T_cell_to_array(cTest2);
    clear cTest2
    final_T = heuristic_hosvd_imputation(T2, cameras_gt2, t_scale);
    res = residual_hosvd(final_T, cameras_gt2, t_scale);
    disp('HO-rSTE');
    disp('Relative Cam Difference | Mean Location Error | Median Location Error | Mean Rotation Error | Median Rotation Error |');
    fprintf('%.4f %.4f %.4f %.4f %.4f \n', res.res2, res.errT_mean, res.errT_median, res.errR_mean, res.errR_median);

else
    final_T = heuristic_hosvd_imputation(T, cameras_gt, t_scale);
    res = residual_hosvd(final_T, cameras_gt, t_scale);
    disp('HO-rSTE');
    disp('Relative Cam Difference | Mean Location Error | Median Location Error | Mean Rotation Error | Median Rotation Error |');
    fprintf('%.4f %.4f %.4f %.4f %.4f \n', res.res2, res.errT_mean, res.errT_median, res.errR_mean, res.errR_median);
end
% save("EPFL_code/"+dataset+"_final_tensor.mat", "T", "cameras_gt", "final_T")

% load("EPFL_code/"+dataset+"_final_tensor.mat")



% save(dataset+"_cameras_est.mat")



function res = residual_hosvd(cur_t, np, t_scale)
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
end

function npmat = normalizenpmat4(npmat)
    n = size(npmat,1)/3;
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
