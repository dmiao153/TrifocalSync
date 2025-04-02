function [errR_mean, rotangles] = new_CompareRotations(camestR,camorigR)

n = size(camestR,3);
rotangles = zeros(1,n);
for i = 1:n
    curest = camestR(:,:,i);
    curorig = camorigR(:,:,i);
    anglemat = curest * curorig';
    anglerad = acos(max(min((trace(anglemat)-1)/2, 1), -1));
    rotangles(i) = anglerad/pi*180;
end

errR_mean = mean(rotangles);
% s = eig(Rb * Rbgt');
% degree1 = acos(real(s(1)));
% s = eig(Rc * Rcgt');
% degree2 = acos(real(s(1)));
% errR_mean = mean([0,abs(degree1),abs(degree2)]);


end