function [cs, cs2] = retrieve_cameras_hosvd(T)
[U,S,sv] = mlsvd(T);
% [U,S,sv] = regularized_hoste(T,[4,4,6],15);

cs = U{1}(:,1:4);
cs2 = U{2}(:,1:4);
% C = U{3}(:,1:6);
% G = S(1:4, 1:4, 1:6);

end

