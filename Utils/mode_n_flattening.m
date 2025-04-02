function flat = mode_n_flattening(T, n)
    ts = size(T);
    % mode1
    if n == 1  
        flat = reshape(T,[ts(1),prod(ts(1:end ~= 1))]);
    end
    
    % mode 2
    if n==2
%         nt = sym(zeros(ts(2),ts(1),ts(3)));
%        expression nt(ts(2),ts(1),ts(3));
        nt = zeros(ts(2),ts(1),ts(3));

        for i = 1:ts(1)
            for j = 1:ts(2)
                for k = 1:ts(3)
                    nt(j,i,k) = T(i,j,k);
                end
            end
        end
        flat = reshape(nt, [ts(2),prod(ts(1:end ~= 2))]);
    end
    
    % mode 3
    if n==3
        % nt = sym(zeros(ts(3),ts(1),ts(2)));

        nt = zeros(ts(3),ts(1),ts(2));

        for i = 1:ts(1)
            for j = 1:ts(2)
                for k = 1:ts(3)
                    nt(k,i,j) = T(i,j,k);
                end
            end
        end
        flat = reshape(nt, [ts(3),prod(ts(1:end ~= 3))]); 
    end
    
end
