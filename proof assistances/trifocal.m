function [T] = trifocal(A,B,C)
% Computes the trifocal tensor for cameras A,B,C 
%   using procedure in http://vigir.missouri.edu/~gdesouza/Research/Conference_CDs/IEEE_CVPR_2009/data/papers/0161.pdf



[L,S,H] = svd(A);
S=S(1:3,1:3);
Aprime = [eye(3) zeros(3,1)];
Bprime = B*H;
Cprime = C*H;
Tprime = zeros(3,3,3);
for i=1:3
    Tprime(i,:,:) = Bprime(:,i)*transpose(Cprime(:,4)) - Bprime(:,4)*transpose(Cprime(:,i));
end
D = (S^(-1))*transpose(L);
T = zeros(3,3,3);
for j=1:3
    T(j,:,:) = D(1,j)*Tprime(1,:,:) + D(2,j)*Tprime(2,:,:) + D(3,j)*Tprime(3,:,:);
end
end

