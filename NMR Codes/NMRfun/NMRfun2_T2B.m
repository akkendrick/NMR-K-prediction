function [loglike, lkpred, likelihood] = NMRfun2(x, Dk, phi, T2, m, n)
% this function computes the likelihood for Kpredicted given K_measure
        logb = x(1); 
        sig = x(2); 

        T2B = x(3);
        
        kk = log10(Dk); 
        lphi = log10(phi);
        
        %lkpred = logb + n*log10(T2);                    % without T2B
        %lkpred = logb - n*log10((T2.^-1)-(T2B.^-1));
        %lkpred = logb + m*lphi + n*log10(T2);
        lkpred = logb + m*lphi - n*log10((T2.^-1)-(T2B.^-1));
        
                
        r = kk - lkpred;

        tau = (1./sig.^2); 
        n = length(r); 
        liketerm1 = (-.5*tau)*(r'*r); 
        liketerm2 = (n/2)*log(tau) - (n/2)*log(2*pi);
        loglike = liketerm2 + liketerm1; 
        likelihood = exp(loglike);
end
