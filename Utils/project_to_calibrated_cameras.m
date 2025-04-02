function [R,t,reflecttrue] = project_to_calibrated_cameras(cam)
% The new camera will be R * [I | -t];

    [U,S,V] = svd(cam(1:3,1:3));
    R = U * V';
    reflecttrue = false;
    if det(R) < 0
        reflecttrue = true;
%         R = R * det(R);
        U(:,3) = -U(:,3);
        R = U * V';
%         R = R * diag([1,-1,-1]);
    end
       
    transmap = R * inv(cam(1:3,1:3));
    
    ccam = transmap * cam;
    
    R = ccam(1:3,1:3);
    t = -inv(R) * ccam(1:3,4);
    
    
end

