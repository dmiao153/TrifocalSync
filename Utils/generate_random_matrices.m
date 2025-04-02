%%% generate_random_camera_matrices using rand(3,4)
%%% input: num_of_P, int. Output: np, cell array. 


function np = generate_random_matrices(num_of_P)

np = {};
for i = 1:num_of_P
    rm = zeros(3,4);
    while rank(rm) < 3
        rm = rand(3,4);
    end
    np{i} =rm;
end

end