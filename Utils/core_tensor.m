function cT = core_tensor()
    cT = zeros(4,4,6);

%     One convention for the exterior square
    cT(3,4,1) = 1;
    cT(4,3,1) = -1;

    cT(2,4,2) = -1;
    cT(4,2,2) = 1;

    cT(2,3,3) = 1;
    cT(3,2,3) = -1;

    cT(1,4,4) = 1;
    cT(4,1,4) = -1;

    cT(1,3,5) = -1;
    cT(3,1,5) = 1;

    cT(1,2,6) = 1;
    cT(2,1,6) = -1;

end
