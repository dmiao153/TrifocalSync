%% This function compares two sets of projection matrices by first 
%% normalizing them and bringing them to the same projection frame

% input: np, np2 are two sets of cameras
%        uncali: is true if they are uncalibrated cameras, and false if
%        calibrated cameras
% output: dd is the difference in norms
%         nn is the norm of the normalized camera matrices
%         Hf is the homography found
%         npmat1 is the 3n x 4 matrix after normalization

function [dd,nn,Hf, npmat1, npmat2] = compare_projection_matrices(np,np2, uncali)

    n = max(size(np));
    
%     final_Hf = zeros(4,4);
%     for i = 1:n
%         for j = 1:n
%             if i ~= j
%                 [Hf, Psa] = findHomography({np{i},np{j}}, {np2{i},np2{j}});
%                 [Hf, Psa] = findHomography({np{1},np{2}}, {np2{1},np2{2}});
%                 final_Hf = final_Hf + Hf;
%             end
%         end
%     end
%     
%     Hf = final_Hf / (n^2-n);
%    [Hf, Psa] = findHomography({np{1},np{2}}, {np2{1},np2{2}});
    [Hf, Psa, vscale] = findHomography_large(np, np2);
    
    npmat1 = reorder_np_cell_to_array(np);
    npmat2 = reorder_np_cell_to_array(np2);
    
    for i = 1:n 
        npmat1(3*(i-1)+1:3*i,:) = npmat1(3*(i-1)+1:3*i,:) * Hf / vscale(i);
    end
    
    if ~uncali
        npmat1 = project_to_calibrated(npmat1);
        npmat2 = project_to_calibrated(npmat2);
    end
        
    npmat1 = normalizenpmat(npmat1);
    npmat2 = normalizenpmat(npmat2);
   
    for i = 1:n
        cur_mat1 = npmat1(3*(i-1)+1:3*i, :);
        cur_mat2 = npmat2(3*(i-1)+1:3*i, :);
        if norm(cur_mat1 + cur_mat2) < norm(cur_mat1 - cur_mat2)
            npmat2(3*(i-1)+1:3*i, :) = cur_mat2 * (-1);
        end
    end
    
    dd1 = norm(npmat1-npmat2,'fro');
%     dd2 = norm((-1) * npmat1 - npmat2, 'fro');
    
    nn = norm(npmat1,'fro');
%     dd = min(dd1,dd2);
    dd = dd1;
    
end



function npmat = normalizenpmat(npmat)

    n = size(npmat,1)/3;
    for i = 1:n
       npmat(3*(i-1)+1:3*i,:) = npmat(3*(i-1)+1:3*i,:) / norm(npmat(3*(i-1)+1:3*i,:),'fro'); 
    end

end
