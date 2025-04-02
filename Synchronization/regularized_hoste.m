% n = 5;
% np = generate_random_matrices(n);
% T = generate_block_trifocal_tensor(np, n);
% [U, S, sv] = hoste2(T , [4,4,6]);
% frob(T)
% frob(T-tmprod(S,{U{1}, U{2}, U{3}}, 1:3))

function [U, S, sv] = regularized_hoste(T , core_size, alpha)

% fprintf("rank of mode 2 flattening is %d \n ",rank(tens2mat(T,1)))

n = size(T);
if n <= 150

    [cov, ind]=regularized_STE_selec_lam(tens2mat(T,1),core_size(1), alpha);
    [U,D,V] = svd(cov);
    v1 = V(:,1:core_size(1));
    d1 = diag(D);

    % fprintf("rank of mode 2 flattening is %d \n ", rank(tens2mat(T,2)))

    [cov, ind]=regularized_STE_selec_lam(tens2mat(T ,2),core_size(2), alpha);
    [U,D,V] = svd(cov);
    v2 = V(:,1:core_size(2));
    d2 = diag(D);

    % fprintf("rank of mode 2 flattening is %d \n ",rank(tens2mat(T,3)))

    [cov, ind]=regularized_STE_selec_lam(tens2mat(T,3),core_size(3), alpha);
    [U,D,V] = svd(cov);
    v3 = V(:,1:core_size(3));
    d3 = diag(D);

    core = tmprod(T, {v1', v2', v3'}, 1:3);

    U = {v1,v2,v3};
    S = core;
    sv = {d1(1:core_size(1)),d2(1:core_size(2)),d3(1:core_size(3))};
else
    U = cell(1,3);
    sv = cell(1,3);
    parfor i = 1:3
        switch i
            case 1
                weights_cur = weights1;
            case 2
                weights_cur = weights2;
            case 3
                weights_cur = weights3;
            otherwise
                n1 = size(weights1,1);
                weights_cur = ones(n1,n1,n1);
        end
        
        [cov, ind]=regularized_STE_selec_lam(tens2mat(T,i),core_size(i), alpha);
        [U,D,V] = svd(cov);
        v1 = V(:,1:core_size(i));
        U{i} = v1;
        d1 = diag(D);   
        sv{i} = d1(1:core_size(i));
    end
        
    S = tmprod(T, {U{1}', U{2}', U{3}'}, 1:3);
    
end
end




