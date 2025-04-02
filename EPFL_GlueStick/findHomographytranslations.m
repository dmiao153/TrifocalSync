% This function maps Ps to Pso through Hf, so Ps{1}* Hf = Pso{1}, Ps{2} *
% Hf = Pso{2}. 
function [ Hf ,Psa, vscales] = findHomographytranslations( ts,tso )
n = size(ts,2);
normPs = cell(n,1);
normPso = cell(n,1);
for i = 1:n
    normPs{i} = ts(:,i) / norm(ts(:,i));
    normPso{i} = tso(:,i) / norm(tso(:,i));
end


csmatrix=zeros(4*n,16+n);

for i = 1:n
    csmatrix(4*(i-1)+1:4*i, 1:16) = rearrange_block_row(normPs{i});
    csmatrix(4*(i-1)+1:4*i, 16+i) = normPso{i}(:);
end

cs=csmatrix;


[u,d,v]=svd(cs);
size(v);
Hf=[v(1,end) v(2,end) v(3,end) v(4,end);v(5,end) v(6,end) v(7,end) v(8,end);
    v(9,end) v(10,end) v(11,end) v(12,end);v(13,end) v(14,end) v(15,end) v(16,end)];
Psa=zeros(n,4);
for i=1:n
    Psa(i,:)=ts(:,i)'*Hf;
%     Pso{i}*v(16+i,end);
end

vscales = v(17:end,end);
% hfp = Hf';
% cs * [hfp(:);v(17,end);v(18,end)];

% a = v(17,end);
% b = v(18,end);

% v = lsqr(cs, zeros(18,1));
% Hf=[v(1,end) v(2,end) v(3,end) v(4,end);v(5,end) v(6,end) v(7,end) v(8,end);
%     v(9,end) v(10,end) v(11,end) v(12,end);v(13,end) v(14,end) v(15,end) v(16,end)];
% Psa=cell(length(Ps),1);
% for i=1:length(Ps)
%     Psa{i}=Ps{i}*Hf;
%     Pso{i}*v(16+i,end);
% end

end
% function  H = findHomography(pp, cc)
%     p1 = pp{1} / norm(pp{1},'fro');
%     p2 = pp{2} / norm(pp{2},'fro');
%     
%     c1 = cc{1} / norm(cc{1},'fro');
%     c2 = cc{2} / norm(cc{2},'fro');
%     
%     P = [p1; p2];
%     C = [c1; c2];
% 
%     
%     reshapedP = kron(eye(4), P);
%     reshapedC = C(:);
% 
%     rref([reshapedP, reshapedC])
%     H = lsqr(reshapedP, reshapedC, 1e-14);
%     H = reshape(H, [4,4]);
%     
% end

% from GpSfM

function csnosym = rearrange_block_row(Ps_s1)
csnosym = zeros(4, 16);
csnosym(1,1)=-Ps_s1(1);
csnosym(2,2)=-Ps_s1(1);
csnosym(3,3)=-Ps_s1(1);
csnosym(4,4)=-Ps_s1(1);


csnosym(1,5)=-Ps_s1(2);
csnosym(2,6)=-Ps_s1(2);
csnosym(3,7)=-Ps_s1(2);
csnosym(4,8)=-Ps_s1(2);


csnosym(1,9)=-Ps_s1(3);
csnosym(2,10)=-Ps_s1(3);
csnosym(3,11)=-Ps_s1(3);
csnosym(4,12)=-Ps_s1(3);


csnosym(1,13)=-Ps_s1(4);
csnosym(2,14)=-Ps_s1(4);
csnosym(3,15)=-Ps_s1(4);
csnosym(4,16)=-Ps_s1(4);

end


function H = project_Hf_to_SO3(Hf)

    H = zeros(4,4);
    H(4,4) = Hf(4,4);
    
    [U,S,V] = svd(Hf(1:3,1:3));
    H(1:3,1:3) = U*V';
    H(1:3,4) = inv(S) * Hf(1:3,4);

end










