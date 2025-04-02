function stacked_exterior_square = calculate_stacked_exterior_square_from_mat(B,n)

    stacked_exterior_square = zeros(3*n, 6);
    
    for i=1:n
        stacked_exterior_square(3*(i-1)+1:3*i, :) = exterior_square(B(3*(i-1)+1:3*i, :));
    end
    
end