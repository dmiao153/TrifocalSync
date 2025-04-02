%% Extracting Epipoles from a single trifocal tensor, using method in Zisserman. 
function [e,ep] = e_from_T(t)

    u = zeros(3,3);
    v = zeros(3,3);

    for i = 1:3

        [a,~,c] = svd(t(:,:,i));
        u(:,i) = a(:,3);
        v(:,i) = c(:,3);

    end



    [u1,d1,v1] = svd(u');
%     d1 = diag(d1);
%     d1(3) = 0;
%     d1 = diag(d1);
%     r2ut = u1 * d1 * v1';
    
    [u2,d2,v2] = svd(v');
%     d2 = diag(d2);
%     d2(3) = 0;
%     d2 = diag(d2);
%     r2vt = u2 * d2 * v2';
    
    
%     e = null(r2ut);
%     ep = null(r2vt);
    e = v1(:,3);
    ep = v2(:,3);

    
end