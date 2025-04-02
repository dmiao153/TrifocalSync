%%% new_T_from_P using the normalization method described in Hartley and
%%% Zisserman in Chapter 15
%%  np: a cell array consisting of the projection matrices. 
%%      there should only be three projection matrices. 

function T = new_T_from_P(np)
    P1 = np{1};
    P2 = np{2};
    P3 = np{3};

    [U,S,V] = svd(P1);    
    S1 = S(:,1:3);

    D = inv(S1)*U';
    np1 = inv(S1) * U' * P1 * V;
    np2 = P2 * V;
    np3 = P3 * V;

    T = zeros(3,3,3);
    % D1 = zeros(3,3,3);
    for i=1:3
        T(:,:,i) = np2(:,i) * np3(:,4)' - np2(:,4) * np3(:,i)';
        % D1(i,:,:) = D;
    end

    T_final = zeros(3,3,3);
    for i=1:3
        for j=1:3
            T_final(:,:,i) = T_final(:,:,i) + T(:,:,j) * D(j,i);
        end
    end

    T = T_final;
end



%% code archive
 
%     rank(T(:,:,1))
%     rank(T(:,:,2))
%     rank(T(:,:,3))
%     rank(reshape(T(1,:,:),[3 3]))
%     rank(reshape(T(2,:,:),[3 3]))
%     rank(reshape(T(3,:,:),[3 3]))
% rng(3)
% s = rng;
% n=3;
% 
% % calculating random rotations
% q = randrot(n,1);
% 
% rot_matrices = {};
% for i=1:n
%     cur_rot = quat2rotm(q(i));
%     rot_matrices = [rot_matrices, cur_rot];
% end
% 
% 
% % constructing the calibrated camera matrices
% % stored in a cell array with n entries: np
% 
% t = rotatepoint(q,[0, 0, 100]);
% mu = [0 0 0];
% Sigma = [20 20 20];
% R = mvnrnd(mu,Sigma,n);
% t = t + R;
% 
% 
% np = {};
% 
% % K = [4.5539e+03,0,1.4164e+02; 0,4.6161e+03,1.1082e+02; 0,0,1];
% for i=1:n
%     T = [eye(3), -transpose(t(i,:))];
%     cur_c = K * rot_matrices{i} * T;
%     np = [np, cur_c];
% end




