%% Retrieve 3 cameras from a trifocal tensor, using the method in Zisserman. 
function [P1,P2,P3] = P_from_T(t)
    P1 = [eye(3),zeros(3,1)];
    
    [e,ep] = e_from_T(t);
    
    e = e/norm(e);
    ep = ep/norm(ep);
    
    P2 = zeros(3,4);
    P3 = zeros(3,4);
    
    p3nmlzr = ep * ep' - eye(3);
    for i = 1:3
        P2(:,i) = t(:,:,i) * ep;
        P3(:,i) = p3nmlzr*(transpose(t(:,:,i)))*e; 
    end
    
    P2(:,4) = e;
    P3(:,4) = ep;
    

end