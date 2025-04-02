% adopted from A New Rank Constraint on Multi-view Fundamental Matrices, and its Application to Camera Location Recovery 
% on the website   https://www.cs.unc.edu/~ronisen/
function [errR_mean, errR_median, errT_mean, errT_median, SO3Mats_corr2, camorigR, t_fit2, torig_centered] = comparison_for_calibrated_cameras(camest, camorig)

n = max(size(camest));
camestR = zeros(3,3,n);
camorigR = zeros(3,3,n);
camestt = zeros(3,n);
camorigt = zeros(3,n);

[Hf, Psa, vscales] = findHomography_large(camest, camorig);
for i = 1:n 
    camest{i} = camest{i} * Hf / vscales(i);
end

for i=1:n
%     [Kest, Rest, test] = vgg_KR_from_P(camest{i});
%     [Korig, Rorig, torig] = vgg_KR_from_P(camorig{i});
    [Rest, test, reflecttrue] = project_to_calibrated_cameras(camest{i});
%     if ismember(-1, diag(Korig))
%         Kest
%         Korig
%         fprintf("buggy camera calibration\n");
%     end


    Rorig = camorig{i}(1:3,1:3) / norm(camorig{i}(1:3,1:3));
    torig = -inv(Rorig) * camorig{i}(1:3,4);
    
%     if det(Rorig) < 0 || reflecttrue
%         Rest
%         Rorig
%         disp("check")
%     end
    
    
    camestR(:,:,i) = Rest;
    camorigR(:,:,i) = Rorig;
    camestt(:,i) = test;
    camorigt(:,i) = torig;
end

[SO3Mats_corr2, MSE_rots2, R_global] = GlobalSOdCorrectLeft(camestR, camorigR);

torig_centered = camorigt - mean(camorigt, 2);
camestt = camestt - mean(camestt,2);

[t_fit2,t_opt2,c_opt2,NRMSE_LUD2,~] = SimpleTransScaleRemove(R_global*camestt,torig_centered,'L1');



%% Translation Error

errT_mean=mean(sqrt(sum((t_fit2 - torig_centered).^2)));

errT_median=median(sqrt(sum((t_fit2 - torig_centered).^2)));

 
%% Rotation Error in degrees

[E_res,e,normE]=CompareRotations(SO3Mats_corr2,camorigR);

errR_mean=E_res(1); errR_median=E_res(2);

errRnorm_mean=normE(1); errRnorm_median=normE(2);

 

end


function new_Rest = correct_signs(Rest, Rorig)
    sign_diff = [1,-1,-1;-1,1,-1;-1,-1,1];
    
    cur_signs = sign(Rest .* Rorig);
    
    
end

