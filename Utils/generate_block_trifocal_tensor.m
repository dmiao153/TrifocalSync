%%% generates block trifocal tensor. input: np, cell array. n, int. output,
%%% T, 3n x 3n x 3n tensor. 
function T = generate_block_trifocal_tensor(np,n)

    T = zeros(3*n,3*n,3*n);
    for i=1:n
        for j=1:n
            for k=1:n
                temp_P = {np{k}, np{i}, np{j}};
                % stacking the big tensor
                %% Using determinant calculation          
                T(3*(i-1)+1:(3*i), 3*(j-1)+1:3*j, 3*(k-1)+1:3*k) = T_from_P(temp_P); 
                %% original_VGG implementation            
                % T(3*(i-1)+1:(3*i), 3*(j-1)+1:3*j, 3*(k-1)+1:3*k) = vgg_T_from_P(temp_P);
                %% Calibrating first camera implementation
                % T(3*(i-1)+1:(3*i), 3*(j-1)+1:3*j, 3*(k-1)+1:3*k) = new_T_from_P(temp_P);  
            end
        end
    end


end