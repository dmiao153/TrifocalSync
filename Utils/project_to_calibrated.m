%ucs is a 3nx4 matrix of uncalibrated cameras
function ccs = project_to_calibrated(ucs)

    n = max(size(ucs)) / 3;
    ccs = zeros(size(ucs));
    
    for i = 1:n
        scs = ucs(3*(i-1)+1:3*i, :);
%         [ki,ri,ti] = vgg_KR_from_P(scs);
        [ri, ti] = project_to_calibrated_cameras(scs);
%         ki = diag(ki);
%         ki = diag(sign(ki));
        ccs(3*(i-1)+1:3*i, :) = ri*[eye(3), -ti];
    end
    

end