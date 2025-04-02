function [lambda,mu] = findRankOne(H,K)
%INPUT: H and K are 4 x 4 symmetric matrices
%OUTPUT: lambda and mu are scalars such that lambda*H + mu*K is ``as rank-1 
%as possible"
%NOTE: This is a helper function for calibrate, used to help find the best 
%homography
Coeff = [];
for i1=1:3
    for i2=i1+1:4
        for j1=1:3
            for j2=j1+1:4
                a = det(H([i1,i2],[j1,j2]));
                b = H(i1,j1)*K(i2,j2)+K(i1,j1)*H(i2,j2)-H(i1,j2)*K(i2,j1)-K(i1,j2)*H(i2,j1);
                c = det(K([i1,i2],[j1,j2]));
                Coeff = [Coeff; a b c];
            end
        end
    end
end
[~,~,V]=svd(Coeff);
v = V(:,end);
m = [v(1) v(2); v(2) v(3)];
[U, D] = eig(m);
[~,ind]=max(abs(diag(D)));
u = U(:,ind);
lambda = u(1);
mu = u(2);
end

