function fm = reorder_fundamental_from_iij(Tiij)
%reorders the Tiij block to get fm
    fm = zeros(3,3);
    fm(1,1) = Tiij(2,3,1);
    fm(1,2) = Tiij(3,1,1);
    fm(1,3) = Tiij(1,2,1);
    
      
    fm(2,1) = Tiij(2,3,2);
    fm(2,2) = Tiij(3,1,2);
    fm(2,3) = Tiij(1,2,2);
    
    
    fm(3,1) = Tiij(2,3,3);
    fm(3,2) = Tiij(3,1,3);
    fm(3,3) = Tiij(1,2,3);
    
    
end

