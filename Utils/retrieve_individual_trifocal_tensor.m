function one_tft = retrieve_individual_trifocal_tensor(T,i,j,k)

    one_tft = T(3*(i-1)+1:3*i, 3*(j-1)+1:3*j, 3*(k-1)+1:3*k);
    
end
