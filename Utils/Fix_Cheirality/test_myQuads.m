B=rand(3,4);
M=myQuads(B);
l = rand(5,1);
g = [l(1) 0 0 0;
     0 l(1) 0 0;
     0 0 l(1) 0;
     l(2) l(3) l(4) l(5)];
mat1 = (transpose(B*g))*(B*g);
mat1 = mat1(1:3,1:3);
mat2 = zeros(3,3);
for i=1:3
    for j=1:3
        mat2(i,j) = (transpose(l(1:4)))*(M{i,j})*(l(1:4));
    end
end
mat1;
mat2;
norm(mat1 - mat2) %should be 0... and it is :)


isSym=zeros(3,3);
for i = 1:3
    for j = 1:3
        M{i,j}-transpose(M{i,j});
    end
end
isSym %should be 0 matrix... and it is :)