function [ psi ] = makePotential(w, qu, qv)

    alpha = exp(w) - 1;    
    
    K = length(qu);
    L = length(qv);
    
    psi = zeros(K,L);    
    
    for k = 1:K
        for l = 1:L
            [m, xi] = marginalize(alpha, qu(k), qv(l));
            assert(all(imag(m) == 0));
            
            ent = zeros(4,1);
                        
            totalNeg = sum(m(m < 0));
            assert(abs(totalNeg) < 1e-5, 'Too much negativity!');            
            
            ent(m <= 0) = 0;
            nz = m > 0;
            ent(nz) = -m(nz) .* log(m(nz));
            
            psi(k,l) = -w * xi - sum(ent);
        end
    end
end
