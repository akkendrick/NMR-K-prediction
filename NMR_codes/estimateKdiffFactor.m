function kDiffFactor = estimateKdiffFactor(K,k_estimates)

[row, col] = size(k_estimates);
kDiffFactor = zeros(size(k_estimates));

for a = 1:row
    for b = 1:col
        
        factor1 = K(a)/k_estimates(a,b);
        factor2 = k_estimates(a,b)/K(a);

        kDiffFactor(a,b) = max([factor1 factor2]);
    end
end


end

