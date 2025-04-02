function error=ransac_ReprError(ProjM,Corresp,Points3D)
% computing dimensions
N=size(Corresp,2);
M=size(ProjM,2);

% compute 3D triangulation if necessary
if nargin~=3
    Points3D_est=triangulation3D(ProjM,Corresp);
elseif size(Points3D,1)==3
    Points3D_est=[Points3D;ones(1,N)];
elseif size(Points3D,1)==4
    Points3D_est=Points3D;
end

% convert to affine coordinates and adapt Corresp matrix
if size(Corresp,1)==3*M
    Corresp=reshape(Corresp,3,N*M);
    Corresp=Corresp(1:2,:)./repmat(Corresp(3,:),2,1);
else
    Corresp=reshape(Corresp,2,N*M);
end

% reproject points
P=cell2mat(ProjM.');
Corresp_est=reshape(P*Points3D_est,3,M*N);
Corresp_est=Corresp_est(1:2,:)./repmat(Corresp_est(3,:),2,1);

Corresp_est = reshape(Corresp_est, [2*M,N]);
Corresp = reshape(Corresp, [2*M,N]);
% compute RMS of distances for each point
error= sqrt(mean((Corresp_est-Corresp).^2,1));

end