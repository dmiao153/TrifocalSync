%%% generate_random_projection_cameras: 
%%  num_of_P: the number of projection cameras to generate, integer
%%  visualization: true to visualize the camera centers, boolean
%%  fundamental: true to generate random fundamental matrices, boolean

function [U,V,np,F, t, rot_matrices, K_total] = generate_random_projection_cameras(num_of_P, visualization, fundamental)

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
    t = rotatepoint(q,[0, 0, 100]);
    mu = [0 0 0];
    % Sigma = [5 7 10];
    Sigma = [1 2 3];
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

            K = [rand(1)*1500, rand(1)*1000, rand(1); 0, rand(1)*1500, rand(1); 0, 0, 1];
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
    if visualization == true
        pt = t;
        figure
        scatter3(pt(:,1), pt(:,2), pt(:,3))
        hold on
        plot3(0,0,0,'x')
        axis equal
        hold off
    end

end