function T = stacked_fm_from_np(np)
n = size(np,2);
T = zeros(3*n,3*n);
for i=1:n
    for j=1:n
        % temp_P = {np{i}, np{j}};
        
        % stacking n-view fundamental matrix
        
        T(3*(i-1)+1:(3*i), 3*(j-1)+1:3*j) = vgg_F_from_P(np{j}, np{i}); 
        
        % verification code
        % [U1,S1,V1] = svd(vgg_F_from_P(np{i}, np{j}));
        % diag(S1)
        % rank(vgg_F_from_P(np{i}, np{j}));
    end
end




