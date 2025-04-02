load("EPFL_code/"+dataset+"_two_view_data.mat");

n = data.n; MAXIT = 20; id = 0; rand_id = 0; N = n;

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

%% Run Baseline
%[result,data_mod] = apply_baseline_new(data);
if strcmp(dataset, "HerzP25")
    load("EPFL_code/HerzP25_good_index.mat");
    data = transform_data(data,good_indices);
end
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

% data = transform_data(data,data_mod.CompInds);
[var] = initialize_admm_ndset(result,data,true);
var
%Run IRLS based optimization
tic; [result1] = iterative_R3_irls_ndset(var,data_mod,MAXIT); tt=toc;
fprintf('Time takes for Optimization : %.2f sec\n',tt);
%fix rotation from Baseline
result1.R_out=result.R_out;
%Decompose Essentials into rotation and translation (apply baseline)
evalc('[result_admm_bs,data_admm] = apply_baseline_toadmm_noR_noT(data,result1);')
disp('Done baseline on ADMM with random initizalition');

disp('ADMM random initialization');
disp('ID | rand_id | N | Mean Essential Error | Median Essential Error |  Mean Location Error | Median Location Error | Mean Rotation Error | Median Rotation Error |gm');
fprintf('%d %.2f %.2f %.4f %.4f %.4f %.4f %.4f %.4f\n',id,rand_id,N,result_admm_bs.errorF,result_admm_bs.errorF_median,result_admm_bs.Tran_err_mean,result_admm_bs.Tran_err_median,result_admm_bs.Rot_err_mean,result_admm_bs.Rot_err_median);



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

