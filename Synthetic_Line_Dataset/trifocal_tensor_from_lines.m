%% Calculate Trifocal Tensors from Lines _ Linear
function [T,P1,P2,P3]=trifocal_tensor_from_lines(l,lp,lpp)
% Input Arguments:
% l, lp, lpp - 3 vectors 3xN containing image lines for image 1,2,3
%              respectively
% will need at least 9 lines
N = size(l,2);
A = zeros(3*N, 27);

for i = 1:N
    l1 = l(1,i); l2 = l(2,i); l3 = l(3,i);
    lp1 = lp(1,i); lp2 = lp(2,i); lp3 = lp(3,i);
    lpp1 = lpp(1,i); lpp2 = lpp(2,i); lpp3 = lpp(3,i);
    
    A(3*(i-1)+1,:) = [0,0,0,0,0,0,0,0,0,...
        l3*lpp1*lp1, l3*lpp1*lp2, l3*lpp1*lp3, l3*lpp2*lp1, l3*lpp2*lp2, l3*lpp2*lp3, l3*lpp3*lp1, l3*lpp3*lp2, l3*lpp3*lp3,...
        -l2*lpp1*lp1, -l2*lpp1*lp2, -l2*lpp1*lp3, -l2*lpp2*lp1, -l2*lpp2*lp2, -l2*lpp2*lp3, -l2*lpp3*lp1, -l2*lpp3*lp2, -l2*lpp3*lp3];
    A(3*(i-1)+2,:) = [-l3*lpp1*lp1, -l3*lpp1*lp2, -l3*lpp1*lp3, -l3*lpp2*lp1, -l3*lpp2*lp2, -l3*lpp2*lp3, -l3*lpp3*lp1, -l3*lpp3*lp2, -l3*lpp3*lp3,...
        0,0,0,0,0,0,0,0,0,...
         l1*lpp1*lp1, l1*lpp1*lp2, l1*lpp1*lp3, l1*lpp2*lp1, l1*lpp2*lp2, l1*lpp2*lp3, l1*lpp3*lp1, l1*lpp3*lp2, l1*lpp3*lp3];
    A(3*(i-1)+3,:) = [l2*lpp1*lp1, l2*lpp1*lp2, l2*lpp1*lp3, l2*lpp2*lp1, l2*lpp2*lp2, l2*lpp2*lp3, l2*lpp3*lp1, l2*lpp3*lp2, l2*lpp3*lp3,...
        -l1*lpp1*lp1, -l1*lpp1*lp2, -l1*lpp1*lp3, -l1*lpp2*lp1, -l1*lpp2*lp2, -l1*lpp2*lp3, -l1*lpp3*lp1, -l1*lpp3*lp2, -l1*lpp3*lp3,...
        0,0,0,0,0,0,0,0,0];
end

[~,~,V] = svd(A);
t = V(:,end);
T = reshape(t,3,3,3);


% compute valid tensor
% epipoles
[~,~,V]=svd(T(:,:,1)); v1=V(:,end);
[~,~,V]=svd(T(:,:,2)); v2=V(:,end);
[~,~,V]=svd(T(:,:,3)); v3=V(:,end);
[~,~,V]=svd([v1 v2 v3].'); epi31=V(:,end);

[~,~,V]=svd(T(:,:,1).'); v1=V(:,end);
[~,~,V]=svd(T(:,:,2).'); v2=V(:,end);
[~,~,V]=svd(T(:,:,3).'); v3=V(:,end);
[~,~,V]=svd([v1 v2 v3].'); epi21=V(:,end);

% using matrices' parameters
E=[kron(eye(3),kron(epi31,eye(3))),-kron(eye(9),epi21)];
[U,S,V]=svd(E); Up=U(:,1:rank(E)); Vp=V(:,1:rank(E)); Sp=S(1:rank(E),1:rank(E));
[~,~,V]=svd(A*Up); tp=V(:,end);
t=Up*tp;
a=Vp*inv(Sp)*tp;

P1=eye(3,4);
P2=[reshape(a(1:9),3,3), epi21];
P3=[reshape(a(10:18),3,3), epi31];
T=reshape(t,3,3,3);


end