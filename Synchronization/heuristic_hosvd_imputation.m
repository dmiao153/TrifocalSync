function final_T = heuristic_hosvd_imputation(T,np, t_scale)

n = round(size(T,1)/3);
have_gt = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uncali = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 1
    np = cell(int32(n));
    for i = 1:n
       np{i} = zeros(3,4); 
    end
    have_gt = false;
end

%% We find the mask of observed entries
w_mask = zeros(n,n,n);

for i = 1:n
    for j = 1:n
        for k = 1:n 
            if i~=j || j~=k || i~=k 
                if nnz(T(3*(i-1)+1:3*i, 3*(j-1)+1:3*j, 3*(k-1)+1:3*k)) > 0
                    w_mask(i,j,k) = 1;
                end
            else
                w_mask(i,i,i) = 1;
            end         
        end
    end
end

fprintf("The number of cameras used for synchronization is %d \n", n);

nw_mask = ~w_mask;
fprintf("The sampling rate is %f \n", (nnz(w_mask))/(n^3));

known_positions_vec = reshape(w_mask, [n,n,n]);
unknown_positions_vec = reshape(nw_mask, [n,n,n]);

w_mask = repelem(w_mask,3,3,3);
nw_mask = repelem(nw_mask, 3,3,3);


%% start of the hosvd iterative optimization
cur_iter = 1;
MAX_iter = 100; % was 150
residuals = zeros(MAX_iter,1); 
algo_cost_history = zeros(MAX_iter,1); 
meanRerr = zeros(MAX_iter,1); 
medianRerr = zeros(MAX_iter,1); 
meanTerr = zeros(MAX_iter,1); 
medianTerr = zeros(MAX_iter,1); 
scale_var = zeros(MAX_iter,1);

%% imputation
% curnw = rand(3*n,3*n,3*n);
curnw = zeros(3*n,3*n,3*n);
cur_t = T + nw_mask .* curnw;

%% Find the iterative soft thresholders
[U,S,sv] = mlsvd(T);

l1 = sv{1}(max(floor(3*n/3),4));
l2 = sv{2}(max(floor(3*n/3),4)); 
l3 = sv{3}(max(floor(3*n/3),6));

sv{1} = sv{1}(sv{1} - l1 > 0);
sv{2} = sv{2}(sv{2} - l2 > 0);
sv{3} = sv{3}(sv{3} - l3 > 0);

r1 = nnz(sv{1});
r2 = nnz(sv{2});
r3 = nnz(sv{3});

fprintf("The initial ranks are (%d, %d, %d)\n", r1,r2,r3)

%% Calculate the residuals

if have_gt
    res = residual_hosvd(cur_t, np, t_scale);
    residuals(cur_iter) = res.res2;
    meanRerr(cur_iter) = res.errR_mean;
    medianRerr(cur_iter) = res.errR_median;
    meanTerr(cur_iter) = res.errR_mean;
    medianTerr(cur_iter) = res.errT_median;
end

r1_history = zeros(MAX_iter,1);
r2_history = zeros(MAX_iter,1);
r3_history = zeros(MAX_iter,1);

rank_truncated = zeros(3*n,3*n,3*n);

while cur_iter < MAX_iter    
    cur_t(isnan(cur_t)) = 0;

    %% imputation needed in each iteration
    % prev_t = cur_t;
    cur_t = cur_t .* w_mask + nw_mask .* rank_truncated;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% UNCOMMENT IF USE HOSVD-HT, COMMENT IF USE ROBUST VARIANT
%     [U,S,sv] = mlsvd_rsi(cur_t,[3*n,3*n,3*n]);
% 
%     sv{1} = sv{1}(sv{1} - l1 > 0);
%     sv{2} = sv{2}(sv{2} - l2 > 0);
%     sv{3} = sv{3}(sv{3} - l3 > 0);
%     r1 = max(4,nnz(sv{1}));
%     r2 = max(4,nnz(sv{2}));
%     r3 = max(6,nnz(sv{3}));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMMENT IF USE HOSVD-HT, UNCOMMENT IF USE ROBUST VARIANT

    [U,S,sv] = regularized_hoste(cur_t,[4,4,6],15);
    
    algo_cost_history(cur_iter) = calculate_svd_cost(cur_t);
    
    r1_history(cur_iter) = r1;
    r2_history(cur_iter) = r2;
    r3_history(cur_iter) = r3;
    
    rank_truncated = tmprod(S, {U{1},U{2},U{3}}, 1:3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% UNCOMMENT IF USE HOSVD-HT, COMMENT IF USE ROBUST VARIANT
%     algo_cost_history(cur_iter) = calculate_svd_cost(cur_t);
% 
%     r1_history(cur_iter) = r1;
%     r2_history(cur_iter) = r2;
%     r3_history(cur_iter) = r3;
% 
% 
%     cs = U{1}(:,1:r1);
%     B = U{2}(:,1:r2);
%     C = U{3}(:,1:r3);
%     G = S(1:r1, 1:r2, 1:r3);
%     rank_truncated = tmprod(G, {cs,B,C}, 1:3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    cur_lambdas = lambda_from_admm(rank_truncated,cur_t);  % calculate the scales with norm reconstructed_T/ norm prev_t

    temp_lam = cur_lambdas(:);
    temp_lam = temp_lam(temp_lam~=0);
    scale_var(cur_iter) = var(temp_lam);
    if cur_iter > 70 && ((scale_var(cur_iter) / scale_var(cur_iter-1) >= 100) || ...
            (scale_var(cur_iter) / scale_var(cur_iter-2) >= 100))
        fprintf("Stopped at iteration %d\n", cur_iter);
        break;
    end
         
    expanded_cur_lambdas = repelem(cur_lambdas,3,3,3);
    cur_t = expanded_cur_lambdas .* cur_t .* w_mask + rank_truncated .* nw_mask;  % cur_t is now the scaled T's, unknown blocks imputed with rank truncated tensor
    
    cur_t(isnan(cur_t)) = 0;
    if have_gt
        res = residual_hosvd(cur_t, np, t_scale);
        residuals(cur_iter) = res.res2;
        meanRerr(cur_iter) = res.errR_mean;
        medianRerr(cur_iter) = res.errR_median;
        meanTerr(cur_iter) = res.errT_mean;
        medianTerr(cur_iter) = res.errT_median;
    end
    
    % iteration update
    cur_iter = cur_iter +1;
    if mod(cur_iter,10) == 0
        fprintf("Done with synchronization iteration %d \n", cur_iter);
    end

end 
if have_gt
    figure
    subplot(3,2,1);
    plot(algo_cost_history(5:cur_iter-1))
    ylabel('cost')
    xlabel('Iteration')
    title('Algorithm Cost Progression')

    subplot(3,2,2)
    plotresidual = residuals(1:cur_iter-1);
    plot(plotresidual * 100)
    ylabel('Residual (%)')
    xlabel('Iteration')
    title('Residual plot')
    
    subplot(3,2,3)
    plot(medianRerr(1:cur_iter-1))
    ylabel('median R error')
    xlabel('Iteration')
    title('median R error plot')
    
    
    subplot(3,2,4)
    plot(medianTerr(1:cur_iter-1))
    ylabel('median T error')
    xlabel('Iteration')
    title('median T error plot')
    
    subplot(3,2,5)
    plot(meanRerr(1:cur_iter-1))
    ylabel('mean R error')
    xlabel('Iteration')
    title('mean R error plot')
    
    subplot(3,2,6)
    plot(meanTerr(1:cur_iter-1))
    ylabel('mean T error')
    xlabel('Iteration')
    title('mean T error plot')
    
else
    plot(algo_cost_history(1:cur_iter-1))
    ylabel('cost')
    xlabel('Iteration')
    title('Algorithm Cost Progression')

end

% figure;
% plot(scale_var(1:cur_iter-1))

final_T = cur_t;
end


function cost = calculate_svd_cost(cur_t)
    [U,S] = mlsvd_rsi(cur_t,[4,4,6]);
    rank_truncated = tmprod(S, {U{1},U{2},U{3}}, 1:3);
    cost = frob(cur_t - rank_truncated)^2;
end

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

%%% This function normalizes T such that each block has frobenius norm 1. 
function T = normalize_tt(T)

n = size(T,1)/3;
for i = 1:n
    for j = 1:n
        for k = 1:n
            if i~=j || j~=k || i~=k
                temp = T(3*(i-1)+1:3*i, 3*(j-1)+1:3*j, 3*(k-1)+1:3*k);
                T(3*(i-1)+1:3*i, 3*(j-1)+1:3*j, 3*(k-1)+1:3*k) = temp / frob(temp);
            end
        end
    end
end

end

%%% this function calculates the lambdas as the quotient of the frobenius 
%%% norms of each block in t and nt, lambda = nt / t

function lambda = lambda_from_admm(nt,t)
n = size(nt,1)/3;

lambda = zeros(n,n,n);
if n >= 70
    parfor i = 1:n    
        for j = 1:n
            for k = 1:n
                if i~=j || j~=k || i~=k
                    cur_ht = tens2mat(t(3*(i-1)+1:3*i,3*(j-1)+1:3*j,3*(k-1)+1:3*k),1);
                    cur_t = tens2mat(nt(3*(i-1)+1:3*i,3*(j-1)+1:3*j,3*(k-1)+1:3*k),1);
                    lambda(i,j,k) = trace(transpose(cur_t) * cur_ht) /((frob(cur_ht))^2);

    %               lambda(i,j,k) = ((frob(cur_t))^2) / trace(transpose(cur_ht) * cur_t) ;
                end
            end
        end
    end
else
    for i = 1:n    
        for j = 1:n
            for k = 1:n
                if i~=j || j~=k || i~=k
                    cur_ht = tens2mat(t(3*(i-1)+1:3*i,3*(j-1)+1:3*j,3*(k-1)+1:3*k),1);
                    cur_t = tens2mat(nt(3*(i-1)+1:3*i,3*(j-1)+1:3*j,3*(k-1)+1:3*k),1);
                    lambda(i,j,k) = trace(transpose(cur_t) * cur_ht) /((frob(cur_ht))^2);

    %               lambda(i,j,k) = ((frob(cur_t))^2) / trace(transpose(cur_ht) * cur_t) ;
                end
            end
        end
    end
end

end

function T = project_T_back(T)
    % iii block constraints
    n = size(T,1)/3;
    z = zeros(3*n,3*n,3*n);
    for i = 1:n
        z(3*(i-1)+1:3*i,3*(i-1)+1:3*i,3*(i-1)+1:3*i) = ones(3,3,3);
    end
    
    T(logical(z)) = 0;
        
    
    % iji and jii block constraints

    z = zeros(3,3,3);
    z(1,:,1) = ones(1,3);
    z(2,:,2) = ones(1,3);
    z(3,:,3) = ones(1,3);
    zl = ~logical(z);
    
    z2 = zeros(3,3,3);
    z2(:,1,1) = ones(3,1);
    z2(:,2,2) = ones(3,1);
    z2(:,3,3) = ones(3,1);    
    z2l = ~logical(z2);
    
    
%     for i = 1:n
%         for j = 1:n
%              if i~=j
%                  temp = T(3*(i-1)+1:3*i, 3*(j-1)+1:3*j, 3*(i-1)+1:3*i);
%                  temp(zl) = 0;
%                  avg1 = (temp(1,1,1) + temp(2,1,2) + temp(3,1,3))/3;
%                  avg2 = (temp(1,2,1) + temp(2,2,2) + temp(3,2,3))/3;
%                  avg3 = (temp(1,3,1) + temp(2,3,2) + temp(3,3,3))/3;
%                  
%                  temp(1,1,1) = avg1;
%                  temp(2,1,2) = avg1;
%                  temp(3,1,3) = avg1;
%                  
%                  temp(1,2,1) = avg2;
%                  temp(2,2,2) = avg2;
%                  temp(3,2,3) = avg2;
%                  
%                  temp(1,3,1) = avg3;
%                  temp(2,3,2) = avg3;
%                  temp(3,3,3) = avg3;
%                  
%                  T(3*(i-1)+1:3*i, 3*(j-1)+1:3*j, 3*(i-1)+1:3*i) = temp;
%                  
%                  
%                  temp = T(3*(j-1)+1:3*j, 3*(i-1)+1:3*i, 3*(i-1)+1:3*i);
%                  temp(z2l) = 0;
%                  
%                  avg1 = (temp(1,1,1) + temp(1,2,2) + temp(1,3,3))/3;
%                  avg2 = (temp(2,1,1) + temp(2,2,2) + temp(2,3,3))/3;
%                  avg3 = (temp(3,1,1) + temp(3,2,2) + temp(3,3,3))/3;
%                  
%                  temp(1,1,1) = avg1;
%                  temp(1,2,2) = avg1;
%                  temp(1,3,3) = avg1;
%                  
%                  temp(2,1,1) = avg2;
%                  temp(2,2,2) = avg2;
%                  temp(2,3,3) = avg2;
%                  
%                  temp(3,1,1) = avg3;
%                  temp(3,2,2) = avg3;
%                  temp(3,3,3) = avg3;
%                  
%                  T(3*(j-1)+1:3*j, 3*(i-1)+1:3*i, 3*(i-1)+1:3*i) = temp;
% 
%                  
%              end
%         end
%     end
    
%    % skew-symmetry constraint
%     for i = 1:3*n
%         T(:,:,i) = 0.5 * (T(:,:,i) - transpose(T(:,:,i)));
%     end
    
        
end

function T = project_fundamental_matrix_constraints(T)
    % iij block constraints
    n = size(T,1)/3;
    for i = 1:n
        for j = 1:n
            if i~=j
                Tiij = xk.T(3*(i-1)+1:3*i,3*(i-1)+1:3*i,3*(j-1)+1:3*j);
                fm = reorder_fundamental_from_iij(Tiij);
                [u,v,s] = svds(fm);
                v(3,3) = 0;
                fm = u * v * s';
                
                T(3*(i-1)+1:3*i,3*(i-1)+1:3*i,3*(j-1)+1:3*j) = reorder_iij_from_fundamental(fm);
            end
        end
    end    
end

