% Adopted from https://github.com/LauraFJulia/TFT_vs_Fund
function [R_t_2,R_t_3,Reconst,T,iter, distances]=STETFTPoseEstimation(Corresp,CalM)
%  Input arguments:
%  Corresp  - 6xN matrix containing in each column, the 3 projections of
%             the same space point onto the 3 images.
%  CalM     - 9x3 matrix containing the M calibration 3x3 matrices for 
%             each camera concatenated.
%
%  Output arguments: 
%  R_t_2    - 3x4 matrix containing the rotation matrix and translation 
%             vector [R2,t2] for the second camera.
%  R_t_3    - 3x4 matrix containing the rotation matrix and translation 
%             vector [R3,t3] for the third camera.
%  Reconst  - 3xN matrix containing the 3D reconstruction of the
%             correspondences.
%  T        - 3x3x3 array containing the trifocal tensor associated to 
%             this triplet of cameras.
%  iter     - number of iterations needed in estimation (always 0)

% Normalization of the data
[x1,Normal1]=Normalize2Ddata(Corresp(1:2,:));
[x2,Normal2]=Normalize2Ddata(Corresp(3:4,:));
[x3,Normal3]=Normalize2Ddata(Corresp(5:6,:));

% Model to estimate T: linear equations
% [cov, T] = STE_prepartion_tft(x1,x2,x3);
[cov, T, distances] = STE_prepartion_tft_full(x1,x2,x3);

% linT=linearTFT(x1,x2,x3);

% tensor denormalization
T=transform_TFT(T,Normal1,Normal2,Normal3,1);

% Find orientation using calibration and TFT
[R_t_2,R_t_3]=R_t_from_TFT(T,CalM,Corresp);

% Find 3D points by triangulation
% Reconst=triangulation3D({CalM(1:3,:)*eye(3,4),CalM(4:6,:)*R_t_2,CalM(7:9,:)*R_t_3},Corresp);
% Reconst=Reconst(1:3,:)./repmat(Reconst(4,:),3,1);
% 
% iter=0;
Reconst = 0;
iter = 0;
end