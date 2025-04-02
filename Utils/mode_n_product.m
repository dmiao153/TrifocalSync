function tfp = mode_n_product(T, U, n)
    st = size(T);
    su = size(U);
    
    tf = mode_n_flattening(T,n);
    if n==1
        tfp = reshape(U * tf, [su(1),st(2),st(3)]);
    end
    if n==2
        % tfp = sym(zeros(st(1),su(1),st(3)));
        tfp = zeros(st(1),su(1),st(3));

        for i=1:st(3)
            tfp(:,:,i) = (U * tf(:,st(1)*(i-1)+1:st(1)*i))';
        end  
    end
    if n==3
        tfp = reshape((U * tf)', [st(1),st(2),su(1)]);
    end
    
end
