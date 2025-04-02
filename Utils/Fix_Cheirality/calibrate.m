function [Bcal, Ccal, l] = calibrate(B,C)
%INPUT: B, C are 3 x 4 matricees
%OUTPUT: Bcal, Ccal are 3 x 4 calibrated matrices with the right column of 
%Bcal normalized to norm 1 and l is a 5-dimensional vector
%NOTE1: Bcal and Ccal are exactly calibrated but only approximately 
%equivalent to B and C up to a homography preserving [I 0]
%NOTE2: +/- [l(1) 0 0 0; 0 l(1) 0 0; 0 0 l(1) 0; l(2) l(3) l(4) *] (any *)
%is best such homography

[H,K] = extractSubspace(B,C);
[lambda,mu] = findRankOne(H,K);
L = lambda*H + mu*K;
[V,D] = eig(L);
[d,ind]=max(abs(diag(D)));
l=sqrt(d)*V(:,ind);
l = [l; 1/norm(B(:,4),'fro')];
g = [l(1) 0 0 0;
     0 l(1) 0 0;
     0 0 l(1) 0;
     l(2) l(3) l(4) l(5)];
B = B*g;
[U1,~,V1] = svd(B(:,1:3));
B(:,1:3) = U1*transpose(V1);
B = (sign(det(B(:,1:3))))*B;
C = C*g;
[U1,~,V1] = svd(C(:,1:3));
C(:,1:3) = U1*transpose(V1);
C = (sign(det(C(:,1:3))))*C;
Bcal = B;
Ccal = C;
end

