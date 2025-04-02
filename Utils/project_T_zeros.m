function T = project_T_zeros(T,n)
    % iii block constraints
    z = zeros(3*n,3*n,3*n);
    for i = 1:n
        z(3*(i-1)+1:3*i,3*(i-1)+1:3*i,3*(i-1)+1:3*i) = ones(3,3,3);
    end
    
    T(logical(z)) = 0;
        
    
    % iji and jii block constraints

    z = zeros(3,3,3);
    z(1,:,1) = ones(1,3);
    z(2,:,2) = ones(1,3);
    z(3,:,3) = ones(1,3);
    zl = ~logical(z);
    
    z2 = zeros(3,3,3);
    z2(:,1,1) = ones(3,1);
    z2(:,2,2) = ones(3,1);
    z2(:,3,3) = ones(3,1);    
    z2l = ~logical(z2);
    
    
    for i = 1:n
        for j = 1:n
             if i~=j
                 temp = T(3*(i-1)+1:3*i, 3*(j-1)+1:3*j, 3*(i-1)+1:3*i);
                 temp(zl) = 0;
                 avg1 = (temp(1,1,1) + temp(2,1,2) + temp(3,1,3))/3;
                 avg2 = (temp(1,2,1) + temp(2,2,2) + temp(3,2,3))/3;
                 avg3 = (temp(1,3,1) + temp(2,3,2) + temp(3,3,3))/3;
                 
                 temp(1,1,1) = avg1;
                 temp(2,1,2) = avg1;
                 temp(3,1,3) = avg1;
                 
                 temp(1,2,1) = avg2;
                 temp(2,2,2) = avg2;
                 temp(3,2,3) = avg2;
                 
                 temp(1,3,1) = avg3;
                 temp(2,3,2) = avg3;
                 temp(3,3,3) = avg3;
                 
                 T(3*(i-1)+1:3*i, 3*(j-1)+1:3*j, 3*(i-1)+1:3*i) = temp;
                 
                 
                 temp = T(3*(j-1)+1:3*j, 3*(i-1)+1:3*i, 3*(i-1)+1:3*i);
                 temp(z2l) = 0;
                 
                 avg1 = (temp(1,1,1) + temp(1,2,2) + temp(1,3,3))/3;
                 avg2 = (temp(2,1,1) + temp(2,2,2) + temp(2,3,3))/3;
                 avg3 = (temp(3,1,1) + temp(3,2,2) + temp(3,3,3))/3;
                 
                 temp(1,1,1) = avg1;
                 temp(1,2,2) = avg1;
                 temp(1,3,3) = avg1;
                 
                 temp(2,1,1) = avg2;
                 temp(2,2,2) = avg2;
                 temp(2,3,3) = avg2;
                 
                 temp(3,1,1) = avg3;
                 temp(3,2,2) = avg3;
                 temp(3,3,3) = avg3;
                 
                 T(3*(j-1)+1:3*j, 3*(i-1)+1:3*i, 3*(i-1)+1:3*i) = temp;

                 
             end
        end
    end
        
end