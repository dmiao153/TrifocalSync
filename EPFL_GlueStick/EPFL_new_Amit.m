dataset = "HerzP8";

load("EPFL_GlueStick/"+dataset+"_two_view_data.mat");

data.corr = {};
n = data.n; MAXIT = 20; id = 0; rand_id = 0; N = n;

for i = 1:n
    for j = 1:n
        if nnz(data.E_est{i,j}) > 0
            data.E(3*(i-1)+1:3*i, 3*(j-1)+1:3*j) = data.E(3*(i-1)+1:3*i, 3*(j-1)+1:3*j)/ norm(data.E(3*(i-1)+1:3*i, 3*(j-1)+1:3*j));
            data.E_est{i,j} = data.E_est{i,j} / norm(data.E_est{i,j});
        else
            data.E(3*(i-1)+1:3*i, 3*(j-1)+1:3*j) = zeros(3,3);
            data.corr{i,j} = [];
        end
    end
end
data.E_gt = data.E;
% evalc('[result,data_mod] = apply_baseline_new(data);')
% disp('Done Baseline');
% disp('Mean Essential Error | Median Essential Error |  Mean Location Error | Median Location Error | Mean Rotation Error | Median Rotation Error |gm');
% disp('Baseline');
% fprintf(' %.4f %.4f %.4f %.4f %.4f %.4f\n',100*result.errorF,100*result.errorF_median,100*result.Tran_err_mean,100*result.Tran_err_median,result.Rot_err_mean,result.Rot_err_median);
% path_to_data=strcat('EPFL_code/Data/',dataset,'/',dataset,'_dense/');
% [Ks, Rs, ts, ccams, Tgt, imsizes] = EPFL_ground_truth_imsizes(im_names, path_to_data);


%% Run Baseline
%[result,data_mod] = apply_baseline_new(data);
evalc('[result,data_mod] = apply_baseline_new(data);')
disp('Done Baseline');

% result
% data_mod

%% Our Method %%
%initialize
[var] = initialize_admm_ndset(result,data_mod, false);
var
%Run IRLS based optimization
tic; [result1] = iterative_R3_irls_ndset(var,data_mod,MAXIT); tt=toc;
fprintf('Time takes for Optimization : %.2f sec\n',tt);
%fix rotation from Baseline
result1.R_out=result.R_out;
%Decompose Essentials into rotation and translation (apply baseline)
evalc('[result_admm_bs,data_admm] = apply_baseline_toadmm_noR_noT(data_mod,result1);')
disp('Done baseline on ADMM');

%% Result
disp('ID | rand_id | N | Mean Essential Error | Median Essential Error |  Mean Location Error | Median Location Error | Mean Rotation Error | Median Rotation Error |gm');
disp('Baseline');
 fprintf('%d %.2f %.2f %.4f %.4f %.4f %.4f %.4f %.4f\n',id,rand_id,N,result.errorF,result.errorF_median,result.Tran_err_mean,result.Tran_err_median,result.Rot_err_mean,result.Rot_err_median);
disp('ADMM');
fprintf('%d %.2f %.2f %.4f %.4f %.4f %.4f %.4f %.4f\n',id,rand_id,N,result_admm_bs.errorF,result_admm_bs.errorF_median,result_admm_bs.Tran_err_mean,result_admm_bs.Tran_err_median,result_admm_bs.Rot_err_mean,result_admm_bs.Rot_err_median);

%%%%% JUST FOR NON PARALLEL RIGID DATASET
[data, newinds] = transform_data_2view(data,data_mod);

% data = transform_data(data,data_mod.CompInds);
[var] = initialize_admm_ndset(result,data,true);
var
%Run IRLS based optimization
tic; [result1] = iterative_R3_irls_ndset(var,data_mod,MAXIT); tt=toc;
fprintf('Time takes for Optimization : %.2f sec\n',tt);
%fix rotation from Baseline
result1.R_out=result.R_out;
%Decompose Essentials into rotation and translation (apply baseline)
evalc('[result_admm_bs,data_admm,t_orig2_cntrd] = apply_baseline_toadmm_noR_noT(data,result1);')
disp('Done baseline on ADMM with random initizalition');

disp('ADMM random initialization');
disp('ID | rand_id | N | Mean Essential Error | Median Essential Error |  Mean Location Error | Median Location Error | Mean Rotation Error | Median Rotation Error |gm');
fprintf('%d %.2f %.2f %.4f %.4f %.4f %.4f %.4f %.4f\n',id,rand_id,N,result_admm_bs.errorF,result_admm_bs.errorF_median,result_admm_bs.Tran_err_mean,result_admm_bs.Tran_err_median,result_admm_bs.Rot_err_mean,result_admm_bs.Rot_err_median);


%%%%%%%%%%%%% Run BATA here %%%%%%%%%%%%%


evalc('[result_admm_bs,data_admm] = apply_baseline_toadmm_noR_noT_BATA(data_mod,result1);')
disp('Done baseline on BATA');

disp('BATA results');
disp('ID | rand_id | N | Mean Essential Error | Median Essential Error |  Mean Location Error | Median Location Error | Mean Rotation Error | Median Rotation Error |gm');
fprintf('%d %.2f %.2f %.4f %.4f %.4f %.4f %.4f %.4f\n',id,rand_id,N,result_admm_bs.errorF,result_admm_bs.errorF_median,result_admm_bs.Tran_err_mean,result_admm_bs.Tran_err_median,result_admm_bs.Rot_err_mean,result_admm_bs.Rot_err_median);




% Translation Averaging
% param.delta = 10^-6;
% param.numofiterinit = 10;
% param.numofouteriter = 10;
% param.numofinneriter = 10;
% param.robustthre = 10^-1;
% 
% tij_index = [];
% tij_observe = [];
% for i = 1:data.n
%     for j = 1:data.n
%         if i < j && length(data.tijGT{i,j}) > 0
%             tij_index = [tij_index, [i;j]];
%             cur_tdata = (result_admm_bs.TOut(:,i) - result_admm_bs.TOut(:,j)) ...
%                 / norm((result_admm_bs.TOut(:,i) - result_admm_bs.TOut(:,j)));
%             tij_observe = [tij_observe, cur_tdata];
%         end
%     end
% end
% 
% t=BATA(tij_index,tij_observe,param);
% % R = get_t_rotation(t, t_orig2_cntrd);
% 
% t_orig = reshape(data.t, [size(data.t,1), size(data.t,3)]);
% [t_fit2,t_opt2,c_opt2,NRMSE_LUD2,~] = SimpleTransScaleRemove(t, t_orig,'L1');
% 
% 
% % Translation Error
% errT_mean=mean(sqrt(sum((t_fit2 - t_orig).^2)));
% errT_median=median(sqrt(sum((t_fit2 - t_orig).^2)));
% 
% fprintf("errT_mean: %f, errT_median: %f\n", errT_mean, errT_median);












function data = transform_data(data,newinds)
IEVM_cut = zeros(1,data.n);
IEVM_cut(newinds) = 1;
IEVM_cut = logical(IEVM_cut);

data.keep = data.keep(IEVM_cut, IEVM_cut);
data.n = sum(IEVM_cut);
data.AdjMat = data.AdjMat(IEVM_cut,IEVM_cut);
data.Hmat = data.Hmat(logical(kron(IEVM_cut,ones(3,1))),logical(kron(IEVM_cut,ones(3,1))));
data.G_gt = data.G_gt(logical(kron(IEVM_cut,ones(3,1))),logical(kron(IEVM_cut,ones(3,1))));
data.corr = data.corr(IEVM_cut,IEVM_cut);
data.R = data.R(:,:,IEVM_cut);
data.t = data.t(:,IEVM_cut);
data.K=data.K(:,:,IEVM_cut);
data.Focal_gt=data.Focal_gt(:,IEVM_cut);

data.CompInds = data.CompInds(IEVM_cut);
data.E=data.E(logical(kron(IEVM_cut,ones(3,1))),logical(kron(IEVM_cut,ones(3,1))));
data.E_est=data.E_est(IEVM_cut,IEVM_cut);
data.tijGT = data.tijGT(IEVM_cut, IEVM_cut);
data.E_gt = data.E;
end


% function R = get_t_rotation(t, t_gt)
% 
% n = size(t,2);
% % direc_est = cell(n,n);
% % direc_gt = cell(n,n);
% cur_rot = zeros(3,3);
% A = [];
% b = [];
% for i = 1:n    
%     for j = 1:n
%         if i < j
%             cur_direc_est = t(:,i) - t(:,j);
%             cur_direc_gt = t_gt(:,i) - t_gt(:,j);
%             
%             ced = cur_direc_est / norm(cur_direc_est);
%             cgtd = cur_direc_gt / norm(cur_direc_gt);
% %             
% %             direct_est{i,j} = cur_direc_est;
% %             direct_gt{i,j} = cur_direc_gt;
%             A = [A; ced(1), ced(2), ced(3), 0,0,0,0,0,0;...
%                 0,0,0,ced(1),ced(2),ced(3),0,0,0;...
%                 0,0,0,0,0,0,ced(1),ced(2),ced(3)];
%             b = [b;cgtd];
%         end
%     end
% end
% 
% R = A\b;
% R = reshape(R,[3,3]);
% 
% end



% function R = get_t_rotation(t, t_gt)
% 
% n = size(t,2);
% direc_est = cell(n,n);
% direc_gt = cell(n,n);
% cur_rot = zeros(3,3);
% for i = 1:n    
%     for j = 1:n
%         if i < j
%             cur_direc_est = t(:,i) - t(:,j);
%             cur_direc_gt = t_gt(:,i) - t_gt(:,j);
%             
%             cur_direc_est = cur_direc_est / norm(cur_direc_est);
%             cur_direc_gt = cur_direc_gt / norm(cur_direc_gt);
%             
%             curv = cross(cur_direc_est,cur_direc_gt);
%             curc = cur_direc_est' * cur_direc_gt;
%             curs = norm(curv);
%             
%             kmat = vgg_contreps(curv);
%             
%             rotation_matrix = eye(3) + kmat + kmat * kmat * ((1-curc)/(curs^2));
%             
%             cur_rot = cur_rot + rotation_matrix; 
%         end
%     end
% end
% 
% 
% R = cur_rot / (n*(n-1)/2);
% 
% [u,s,v] = svd(R);
% s = eye(3);
% R = u * s * v';
% end
% 

