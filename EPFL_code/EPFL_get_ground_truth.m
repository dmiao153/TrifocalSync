% Adopted from https://github.com/LauraFJulia/TFT_vs_Fund
function [K, R_true, t_true, im_names] = EPFL_get_ground_truth(dataset)
% dataset='Herz_P8';  
% dataset='Herz-Jesu-P8';

%% Some parameters

path_to_data=strcat('EPFL_code/Data/',dataset,'/',dataset,'_dense/');

%% Recover correspondances

dir_data = dir(strcat('EPFL_code/Data/',dataset,'/images/*.png'));
n = length(dir_data);
im_names = cell(1,n);
for i = 1:length(dir_data)
   im_names{i} = dir_data(i).name; 
end
    
    

% get ground truth emat and estimate emat
[Ks, R_true, ts, ccams, Tgt, imsizes, t_true] = EPFL_ground_truth_imsizes(im_names, path_to_data);

K = Ks(:,:,1);
save("EPFL_code/ground_truth_"+dataset+".mat", "K", "R_true", "t_true")
% load("EPFL_code/ground_truth_fountainP11.mat");

end