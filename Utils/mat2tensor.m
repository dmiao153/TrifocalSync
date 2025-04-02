function t = mat2tensor(A, n)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    ds = size(A,1);
    if n == 1
        t = reshape(A,[ds,ds,ds]);
    elseif n == 2
         t = reshape(A,[ds,ds,ds]);
         for i = 1:ds
             t(:,:,i) = transpose(t(:,:,i));
         end
    elseif n==3
        t = zeros(ds,ds,ds);
        for i = 1:ds
            t(:,:,i) = reshape(A(i,:),[ds,ds]);
        end
            
    end
end

