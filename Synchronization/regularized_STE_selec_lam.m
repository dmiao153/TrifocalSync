function [cov, ind]=regularized_STE_selec_lam(X,dd, alpha)
[D,N]=size(X);
initcov=eye(D);
eps=10^-10;
Lam = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15];
% Lam = [2,4,8,12,15];
% Lam = [8];
normalizer = 1/(1+alpha);
normalizer2 = 1 - normalizer;
% Lam = [2];
for i = 1:length(Lam)
    lam = Lam(i);
    oldcov=initcov-1;
    cov=initcov;
    iter=1;

    while norm((oldcov-cov),'fro')>10^-10 & iter<100
        oldcov=cov;
%         temp = (X'/(cov+eps*eye(D))).*X';
%         w = 1./(sum(temp,2)+eps);
%         cov = X * (w .* (X'))/N*D;

        temp = (X'/(cov+eps*eye(D))).*X';
        w = 1./(sum(temp,2)+eps);
        cov = X * (w .* (X'));
        cov = (normalizer * D / N) * cov + normalizer2 * eye(D);
    
        
        [U,S]=svd(cov);
%         [U,S,V] = randpca(cov, size(cov,1));
        S1=diag(real(S));
        S1((dd+1):end)=mean(S1((dd+1):end))/lam;

        cov = U*diag(S1)*U';
        cov = cov/trace(cov);
        iter=iter+1;
    end
    COV{i} = cov;
    C = U(:,end)*U(:,end)'*X;
    dist = ((sum(C.*C,1)).^.5);
    Dis(i,:) = dist;

end

threshold = median(Dis(:));
inliers = sum(Dis<threshold,2);
[B,ind] = max(inliers);
cov = COV{ind};


end