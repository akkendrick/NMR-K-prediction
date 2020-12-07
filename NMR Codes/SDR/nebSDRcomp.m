% Run SDR boostrap on Nebraska Data
clear
load neb_data_3stk.mat

loggingDepth = neb_T2ML_stack3{:,1};

% Convert to seconds
T2ML = neb_T2ML_stack3{:,2}*10^-3;

WBF_topDepth = neb_WBF_K{:,1};
WBF_bottomDepth = neb_WBF_K{:,2};
WBF_K = neb_WBF_K{:,4}; % in ft/day

% Convert K data to m/s
WBF_K = WBF_K * (3.528*10^-6);

%
load neb_data_3stk_phi.mat
phi = neb_data_3stk_phi{:,3};

SDR_K = @(b,m,n,phi,T2ML) (b.*(phi.^m).*(T2ML).^n);
SDR_b = @(K,m,n,phi,T2ML) K./((phi.^m).*(T2ML).^n);

m = 4;
n = 2;

for kk = 1:length(WBF_topDepth)
    currentDepths = loggingDepth(loggingDepth > WBF_topDepth(kk) &...
        loggingDepth < WBF_bottomDepth(kk));
    currentT2ML = T2ML(loggingDepth > WBF_topDepth(kk) &...
        loggingDepth < WBF_bottomDepth(kk));
    currentPhi = phi(loggingDepth > WBF_topDepth(kk) &...
        loggingDepth < WBF_bottomDepth(kk));
    
    currentWBF_K = ones(length(currentDepths),1).*WBF_K(kk);
    
    currentSDRb = SDR_b(currentWBF_K,m,n,0.4,currentT2ML);
    
    %currentSDRK = SDR_K(
    
    % Follow upscaling equation from Dlubac et al. 2013
    p = length(currentDepths);
    
    intervalT2ML = 1/p*sum(currentT2ML);
    intervalb = 1/p*sum(currentSDRb);
    intervalPhi = 1/p*sum(currentPhi);
    
    bProfile(kk) = intervalb;    
    T2MLProfile(kk) = intervalT2ML;
    phiProfile(kk) = intervalPhi;
end

logK = log10(WBF_K);
logT2ML = log10(T2MLProfile)';
logPhi = log10(phiProfile)';

Nboot = 2000;
%logPhi = ones(length(logK),1);

[b_boot, n_boot, m_boot] = bootstrap_fun([logT2ML, logPhi, logK], Nboot, n, m);   % m, n fixed
