function [errR_mean, errT_mean, calesterr, calgterr] = angle_diff(esttft, gttft)
esttft = esttft / frob(esttft);
gttft = gttft / frob(gttft);
[p1,p2,p3] = P_from_T(esttft);
[Bcal, Ccal, l] = calibrate(p2,p3);

[gtp1,gtp2,gtp3] = P_from_T(gttft);
[bcalgt, ccalgt, l] = calibrate(gtp2,gtp3);

% [~, Rb, tb] = vgg_KR_from_P(Bcal);
% [~, Rc, tc] = vgg_KR_from_P(Ccal);
% 
% [~, Rbgt, tbgt] = vgg_KR_from_P(bcalgt);
% [~, Rcgt, tcgt] = vgg_KR_from_P(ccalgt);

after_calT = T_from_P({p1,Bcal, Ccal});
calesterr = scaled_diff_frobenius_norm(after_calT, esttft);

after_calTgt = T_from_P({gtp1,bcalgt, ccalgt});
calgterr = scaled_diff_frobenius_norm(after_calTgt, gttft);


Rb = Bcal(:,1:3);
tb = -inv(Rb) * Bcal(:,4);

Rbgt = bcalgt(:,1:3);
tbgt = -inv(Rbgt) * bcalgt(:,4);

Rc = Ccal(:,1:3);
tc = -inv(Rc) * Ccal(:,4);

Rcgt = ccalgt(:,1:3);
tcgt = -inv(Rcgt) * ccalgt(:,4);



camestR = zeros(3,3,3);
camorigR = zeros(3,3,3);

camestR(:,:,1) = eye(3,3);
camestR(:,:,2) = Rb;
camestR(:,:,3) = Rc;

camorigR(:,:,1) = eye(3,3);
camorigR(:,:,2) = Rbgt;
camorigR(:,:,3) = Rcgt;

camestt = zeros(3,3);
camorigt = zeros(3,3);

camestt(:,1) = p1(:,4);
camestt(:,2) = tb;
camestt(:,3) = tc;

camorigt(:,1) = gtp1(:,4);
camorigt(:,2) = tbgt;
camorigt(:,3) = tcgt;


[SO3Mats_corr2, MSE_rots2, R_global] = GlobalSOdCorrectLeft(camestR, camorigR);

% s = eig(Rb * Rbgt');
% degree1 = acos(real(s(1)));
% s = eig(Rc * Rcgt');
% degree2 = acos(real(s(1)));

% errR_mean = mean([0,abs(degree1),abs(degree2)]);
% [E_res,e,normE]=CompareRotations(SO3Mats_corr2,camorigR);
% errR_mean=E_res(1); errR_median=E_res(2);
% errRnorm_mean=normE(1); errRnorm_median=normE(2);

[errR_mean, rotangles] = new_CompareRotations(SO3Mats_corr2, camorigR);

torig_centered = camorigt - mean(camorigt, 2);
camestt = camestt - mean(camestt,2);
[t_fit2,t_opt2,c_opt2,NRMSE_LUD2,~] = SimpleTransScaleRemove(R_global*camestt,torig_centered,'L2');

errT_mean=mean(sqrt(sum((t_fit2 - torig_centered).^2)));
errT_median=median(sqrt(sum((t_fit2 - torig_centered).^2)));

% if calesterr > 0.5
%     calesterr
%     errT_mean
%     errR_mean
%     calgterr
% end

end
