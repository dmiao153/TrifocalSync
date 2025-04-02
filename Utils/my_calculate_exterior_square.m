function es = my_calculate_exterior_square(x)

    n = size(x,1) / 3;
    
    col_idx = [1,2;1,3;1,4;2,3;2,4;3,4];
    row_idx = [2,3;1,3;1,2];
    
    
    for k = 1:n
        for i=1:3
            for j=1:6
                temp_x = x(3*(k-1)+1:3*k,:);
                temp_mat = temp_x(row_idx(i,:), col_idx(j,:));
                ess(i,j) = temp_mat(1,1) * temp_mat(2,2) - temp_mat(2,1) * temp_mat(1,2);
                % ess(i,j) = temp_mat(1,1) * ((temp_mat(2,2)*temp_mat(3,3)) - (temp_mat(2,3)*temp_mat(3,2))) - temp_mat(1,2) * ((temp_mat(2,1)*temp_mat(3,3)) - (temp_mat(2,3)*temp_mat(3,1))) + temp_mat(1,3) * ((temp_mat(2,1)*temp_mat(3,2)) - (temp_mat(2,2)*temp_mat(3,1)));
            end
        end

        ess(2,:) = ess(2,:) * (-1);
    
        es(3*(k-1)+1:3*k,1:6) = ess;
    end



end