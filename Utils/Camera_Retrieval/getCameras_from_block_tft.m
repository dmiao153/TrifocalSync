% Gets Projective Cameras from the trifocal tensors. Lets say that at least
% two should be the same camera in the indices. 
% 
function [indices, cm, Ps2, check, uPs] = getCameras_from_block_tft(indices, trifocal_tensors, num_cameras)

% indices: [1,2,3;2,3,4;...;n-2,n-1,n]

% NOTE: The ijk block in the block trifocal tensor has index kij. 

% trifocal_tensors: cell array of the trifocal tensors for corresponding
%                   indices
% num_cameras: total number of cameras
% This requires that for each row in the indices, the first three elements
% must have appeared before. 

cm = {};
Ps2 = {};
uPs = {};
check = zeros(1,num_cameras);


for a = 1:size(indices,1)
    if a == 1
        [pi,pj,pk] = P_from_T(trifocal_tensors{a,1}{1}); % The i,j,k block
        cm(a,1:3) = {pi, pj, pk};
        Ps2(indices(a,1)) = {pi};
        Ps2(indices(a,2)) = {pj};
        Ps2(indices(a,3)) = {pk};
        check(indices(a,1)) = 1;
        check(indices(a,2)) = 1;
        check(indices(a,3)) = 1;
    end
    pi = Ps2{indices(a,1)};
    pj = Ps2{indices(a,2)};
    pk = Ps2{indices(a,3)};
    cm(a,1:3) = {pi, pj, pk};

    
    [p2i,p2j,p2k] = P_from_T(trifocal_tensors{a,2}{1}); % The i,j,k block
    cm(a, 4:6) = {p2i, p2j, p2k};

    
    [C,ia,ib] = intersect(indices(a,1:3), indices(a,4:6));
    
    [Hf, ~] = findHomography({cm{a, ib(1) + 3}, cm{a,ib(2)+3}}, {cm{a,ia(1)}, cm{a,ia(2)}});
    % [Hf, Psa] = findHomography({cm{a,ia(1)}, cm{a,ia(2)}}, {cm{a, ib(1) + 3}, cm{a,ib(2)+3}});
    %% Find the homography that maps p2i and p2j to pj and pk

    Hf = double(Hf);
    Hf = Hf/norm(Hf,'fro');

    uPs{a,1} = pi;
    uPs{a,2} = pj;    
    uPs{a,3} = pk;
    
    uPs{a,4} = p2i * (Hf);
    uPs{a,5} = p2j * (Hf);
    uPs{a,6} = p2k * (Hf); 
    
    newp2s = {uPs{a,4},uPs{a,5},uPs{a,6}};
    for b = 4:6
        if ~check(indices(a,b))
            Ps2{indices(a,b)} = uPs{a,b};
            check(indices(a,b)) = 1;
        end
    end
end


end