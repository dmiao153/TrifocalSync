function es = exterior_square(A)

es = zeros(3,6);
col_idx = [1,2;1,3;1,4;2,3;2,4;3,4];
row_idx = [2,3;1,3;1,2];
for i=1:3
    for j=1:6
        es(i,j) = det(A(row_idx(i,:), col_idx(j,:)));
    end
end

es(2,:) = es(2,:) * (-1);

end

