function Tiij = reorder_iij_from_fundamental(fm)
%reorders the Tiij block to get fm
    Tiij = zeros(3,3,3);
    Tiij(2,3,1) = fm(1,1);
    Tiij(3,1,1) = fm(1,2);
    Tiij(1,2,1) = fm(1,3);
    
    Tiij(:,:,1) = Tiij(:,:,1) - (Tiij(:,:,1))';
    
    Tiij(2,3,2) = fm(2,1);
    Tiij(3,1,2) = fm(2,2);
    Tiij(1,2,2) = fm(2,3);
    
    Tiij(:,:,2) = Tiij(:,:,2) - (Tiij(:,:,2))';
    
    Tiij(2,3,3) = fm(3,1);
    Tiij(3,1,3) = fm(3,2);
    Tiij(1,2,3) = fm(3,3);
    
    Tiij(:,:,3) = Tiij(:,:,3) - (Tiij(:,:,3))';
    
end

