function [K,R,t] = getKRt_from_Ps(Ps)
    n = size(Ps,2);
    K = zeros(3,3,n);
    R = zeros(3,3,n);
    t = zeros(3,n);
    
    for i = 1:n
%         Ps{i} = Ps{i} / norm(Ps{i});
        [K(:,:,i), R(:,:,i), t(:,i)] = vgg_KR_from_P(Ps{i});
    end
    

end