function [Ks, Rs, ts, ccams, Tgt, imsizes, t_real] = EPFL_ground_truth_imsizes(im_names, path_to_data)

    n = size(im_names,2);
    Ks = zeros(3,3,n); Rs = zeros(3,3,n); ts = zeros(3,1,n); ccams = zeros(3,4,n); 
    Tgt = zeros(3*n,3*n,3*n); 
    imsizes = zeros(2,1,n);
    t_real = zeros(3,n);
    for i = 1:n
        [K,R,t,imsize]=readCalibrationOrientation_EPFL(path_to_data,im_names{i});
        Ks(:,:,i) = K; Rs(:,:,i) = R; ts(:,:,i) = t; t_real(:,i) = t;
        ccams(:,:,i) = [R, t];
        imsizes(:,i) = imsize;
    end

    for i = 1:n
        for j = 1:n
            for k = 1:n
                Tgt(3*(i-1)+1:3*i,3*(j-1)+1:3*j,3*(k-1)+1:3*k) = T_from_P({ccams(:,:,k), ccams(:,:,i), ccams(:,:,j)});
            end
        end
    end
    
    
end