% This Script generates a scene of lines, estimates trifocal tensors from
% lines, then synchronizes the trifocal tensors. Might need to
% uncomment/comment some code for the synchronization procedure depending
% on which synchronization variant is used. 
clear all
[U,V,np, F, K_total, n, scene_points, image_points, linepoints] = generate_random_scene(1);


T = zeros(3*n,3*n,3*n);
count = 0;
for i = 1:n
    for j = 1:n
        for k = 1:n
            if i < j && j < k 
                [Test,P1,P2,P3]=trifocal_tensor_from_lines(linepoints(:,:,i),linepoints(:,:,j),linepoints(:,:,k));
                Test = T_from_P({P1,P2,P3});
                Tgt = transform_TFT(T_from_P({np{i}, np{j}, np{k}}), K_total{i}, K_total{j},K_total{k},1);

              
                if sum(sign(Test .* Tgt), 'all') < 0
                    P2 = P2 * (-1);
                    count = count + 1;
                end
                
                T(3*(i-1)+1:3*i,3*(j-1)+1:3*j,3*(k-1)+1:3*k) = transform_TFT(T_from_P({P3,P1,P2}),K_total{k}, K_total{i},K_total{j},1);
                T(3*(j-1)+1:3*j,3*(i-1)+1:3*i,3*(k-1)+1:3*k) = transform_TFT(T_from_P({P3,P2,P1}),K_total{k}, K_total{j},K_total{i},1);
                T(3*(k-1)+1:3*k,3*(j-1)+1:3*j,3*(i-1)+1:3*i) = transform_TFT(T_from_P({P1,P3,P2}),K_total{i}, K_total{k},K_total{j},1);
                T(3*(j-1)+1:3*j,3*(k-1)+1:3*k,3*(i-1)+1:3*i) = transform_TFT(Test,K_total{i}, K_total{j},K_total{k},1);
                T(3*(i-1)+1:3*i,3*(k-1)+1:3*k,3*(j-1)+1:3*j) = transform_TFT(T_from_P({P2,P1,P3}),K_total{j}, K_total{i},K_total{k},1);
                T(3*(k-1)+1:3*k,3*(i-1)+1:3*i,3*(j-1)+1:3*j) = transform_TFT(T_from_P({P2,P3,P1}),K_total{j}, K_total{k},K_total{i},1);
            end
        end
    end
end

[metrics, metrics_rotation, metrics_translation] = compare_block_tensor(T, n, np);
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

final_T = heuristic_hosvd_imputation(T, np, 1);

[res, cameras_est] = residual_hosvd(final_T, np);
disp('HO-rSTE');
disp('Relative Cam Difference | Mean Location Error | Median Location Error | Mean Rotation Error | Median Rotation Error |');
fprintf('%.4f %.4f %.4f %.4f %.4f \n', res.res2, res.errT_mean, res.errT_median, res.errR_mean, res.errR_median);

function [res, cameras_est] = residual_hosvd(cur_t, np)
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
res.errT_mean = errT_mean;
res.errT_median = errT_median;
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

function [U,V,np, F, K_total, n, scene_points, image_points, line_points] = generate_random_scene(seed)

rng(seed)
s = rng; n = 20; sNum = 50;
[U,V,np, F, K_total, t] = generate_random_projection_cameras_temp(n, false, false);


Tgt = generate_block_trifocal_tensor(np, n);

scene_points = (rand(4,sNum)-0.5) * 20;

scene_points(4,:) = 1;
image_points = zeros(3,sNum, n);
lnum = cast(sNum/2, "int32");
line_points = zeros(3, lnum, n);
for i = 1:n
    image_points(:,:,i) = np{i} * scene_points;
    for j = 1:lnum
       lineeq = cross(image_points(:,2*j-1,i), image_points(:,2*j,i));
       line_points(:,j,i) = lineeq / norm(lineeq);
    end
end

%%%Gaussian Pixel Wise Noise
noise = normrnd(0,0.0001, 3, sNum/2, n);
fprintf("The noise level is about: %f \n", frob(noise) / frob(line_points));
line_points = line_points + noise;


end


function [U,V,np,F, K_total, t] = generate_random_projection_cameras_temp(num_of_P, visualization, fundamental)

% function np = generate_random_projection_cameras(num_of_P, visualization, fundamental)
%     rng(2)
%     s = rng;
    n=num_of_P;

    % calculating random rotations, q are the random quaternions, and
    % rotmatrices is a cell of 3x3 rotation matrices
    q = randrot(n,1);

    rot_matrices = {};
    for i=1:n
        cur_rot = quat2rotm(q(i));
        rot_matrices = [rot_matrices, cur_rot];
    end
    % constructing the calibrated camera matrices stored in a cell array 
    % with n entries: np. t are the acquired camera centers
    % t = rotatepoint(q,[0, 0, 50000]);
    t = rotatepoint(q,[0, 100, 0]);
    mu = [0 0 0];
%     Sigma = [5 7 10];
    Sigma = [10 20 30];
    R = mvnrnd(mu,Sigma,n);
    t = t + R;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Centering the cameras
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t_avg = sum(t)/n;    
    t = t-t_avg;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    np = {};
    % K = [4.5539e+03,0,1.4164e+02; 0,4.6161e+03,1.1082e+02; 0,0,1];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Verifying low rank Amit Singer F=A+A^T
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    U = zeros(3*n,3);
    V = zeros(3*n,3);
    F = zeros(3*n,3*n);
    K_total = {};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    K = eye(3);
    for i=1:n
        T = [eye(3), -transpose(t(i,:))];
        if fundamental == true
            % K = [rand(1)*2, rand(1)*2, rand(1)*2; 0, rand(1)*2, rand(1)*2; 0, 0, 1];
            K = [rand(1)*150, 0, rand(1)*100; 0, rand(1)*150, rand(1)*100; 0, 0, 1];
        end
        cur_c = K * rot_matrices{i} * T;
        np = [np, cur_c];
        K_total = [K_total, K];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        U(3*(i-1)+1:3*i,:) = inv(K)' * rot_matrices{i} * vgg_contreps(t(i,:));
        V(3*(i-1)+1:3*i,:) = inv(K)' * rot_matrices{i};
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:n
        for j=1:n
            F(3*(i-1)+1:3*i, 3*(j-1)+1:3*j) = inv(K_total{i})'* rot_matrices{i} * (vgg_contreps(t(i,:)) - vgg_contreps(t(j,:))) * rot_matrices{j}'*inv(K_total{j});
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% plotting the camera center positions -- optional
    cam1vec = rot_matrices{1} * [0;100;0];
    points2 = 30 * det(rot_matrices{1}) * rot_matrices{1}(3,1:3);
    if visualization == true
        pt = t;
        figure
        scatter3(pt(:,1), pt(:,2), pt(:,3))
        hold on
        plot3([pt(1,1) + 0,pt(1,1) + points2(1)], [pt(1,2) + 0,pt(1,2) + points2(2)], [pt(1,3) + 0,pt(1,3) + points2(3)])
        plot3(pt(1,1), pt(1,2), pt(1,3), 'x')
        plot3(0,0,0,'x')
        axis equal
        hold off
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
            if nnz(cTest(3*(i-1)+1:3*i, 3*(j-1)+1:3*j, 3*(k-1)+1:3*k)) > 0
                metrics(i,j,k) = scaled_diff_frobenius_norm(cTest(3*(i-1)+1:3*i, 3*(j-1)+1:3*j, 3*(k-1)+1:3*k), T_gt(3*(i-1)+1:3*i, 3*(j-1)+1:3*j, 3*(k-1)+1:3*k));
                [errR_mean, errT_mean] = angle_diff(cTest(3*(i-1)+1:3*i, 3*(j-1)+1:3*j, 3*(k-1)+1:3*k), T_gt(3*(i-1)+1:3*i, 3*(j-1)+1:3*j, 3*(k-1)+1:3*k));
                metrics_rotation(i,j,k) = errR_mean;
                metrics_translation(i,j,k) = errT_mean;
            end
        end
    end
end


end






