%% This code is used to change the final_data cell array tensor to a normal numerical array
function T = final_T_cell_to_array(final_data)
n = size(final_data,1);
T = zeros(3*n,3*n,3*n);
for i = 1:n
    for j = 1:n
        for k = 1:n
            if size(final_data{i,j,k},1)~= 0
                T(3*(i-1)+1:3*i,3*(j-1)+1:3*j,3*(k-1)+1:3*k) = final_data{i,j,k};
            end
        end
    end
end

end