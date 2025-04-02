function [M] = myQuads(B)
%INPUT: B is 3 x 4 matrix 
%OUTPUT: M is 3 x 3 cell array of 4 x 4 symmetric matrices
%NOTE: Expresses entries of left 3 x 3 submatrix of (transpose(B*g))*(B*g)
%where g=[l1 0 0 0; 0 l1 0 0; 0 0 l1 0; l2 l3 l4 l5]
%as quadratic forms in l1, l2, l3, l4

G = transpose(B)*B;
M = cell(3,3);
for i = 1:3
    for j = 1:3
        Q = zeros(4,4);
        Q(1,1) = Q(1,1)+G(i,j);
        Q(1,j+1) = Q(1,j+1)+G(i,4);
        Q(1,i+1) = Q(1,i+1)+G(4,j);
        Q(i+1,j+1) = Q(i+1,j+1)+G(4,4);
        Q = (1/2)*(Q + transpose(Q));
        M{i,j} = Q;
    end
end
end




