function [H,K] = extractSubspace(B,C)
%INPUT B, C are 3 x 4 matrices
%OUTPUT H, K are 4 x 4 symmetric matrices corresponding to minimal right
%singular vectors of linearized quadratic system for l(1:4) parameterizing
%homography

B = B/norm(B,'fro');
C = C/norm(C,'fro');
M = myQuads(B);
N = myQuads(C);

for i=1:3
    for j=1:3
        Q = M{i,j};
        M{i,j} = [Q(1,1) 2*Q(1,2) 2*Q(1,3) 2*Q(1,4) Q(2,2) 2*Q(2,3) 2*Q(2,4) Q(3,3) 2*Q(3,4) Q(4,4)];
        Q = N{i,j};
        N{i,j} = [Q(1,1) 2*Q(1,2) 2*Q(1,3) 2*Q(1,4) Q(2,2) 2*Q(2,3) 2*Q(2,4) Q(3,3) 2*Q(3,4) Q(4,4)];
    end
end

CoeffB = [M{1,1}-(1/3)*(M{1,1}+M{2,2}+M{3,3}); M{1,2}; M{1,3}; M{2,1}; M{2,2}-(1/3)*(M{1,1}+M{2,2}+M{3,3}); M{2,3}; M{3,1}; M{3,2}; M{3,3}-(1/3)*(M{1,1}+M{2,2}+M{3,3})];
CoeffC = [N{1,1}-(1/3)*(N{1,1}+N{2,2}+N{3,3}); N{1,2}; N{1,3}; N{2,1}; N{2,2}-(1/3)*(N{1,1}+N{2,2}+N{3,3}); N{2,3}; N{3,1}; N{3,2}; N{3,3}-(1/3)*(N{1,1}+N{2,2}+N{3,3})];
Coeff = [CoeffB; CoeffC];
[~,~,V] = svd(Coeff);
h = V(:,end-1);
k = V(:,end);
H = [h(1) h(2) h(3) h(4);
     h(2) h(5) h(6) h(7);
     h(3) h(6) h(8) h(9);
     h(4) h(7) h(9) h(10)];
K = [k(1) k(2) k(3) k(4);
     k(2) k(5) k(6) k(7);
     k(3) k(6) k(8) k(9);
     k(4) k(7) k(9) k(10)]; 
end

