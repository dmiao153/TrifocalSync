function r = my_kron(q,w)
    %% q has size n,1
    %% w has size 3,4
    
    nn = size(q,1);
    r = q(1) * w;
    for i = 2:nn
        r = [r;q(i) * w];
    end
    
end