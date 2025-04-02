function metric_res = scaled_diff_frobenius_norm(esttft, gttft)
    esttft = esttft / frob(esttft);
    gttft = gttft / frob(gttft); 
    metric_res = min(frob(esttft-gttft), frob(esttft+gttft));
end
