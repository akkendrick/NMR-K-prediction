% Analyze SOE by plotting K_SOE vs K_DPP
%baseDir = '/Volumes/GoogleDrive/My Drive/Stanford/USGS Project/USGS Data/';
baseDir = 'I:\My Drive\Stanford\USGS Project\Field Data\USGS Data\';

n = [2 1 0];

k_estimates = [];
k_names = {'1:1','SOE n=2','SOE n=1','SOE n=0'};
k_sym = {'+','*','o'};

for k = 1:length(n)
    k
    site = 'Site2-WellPN2'
    siteName = site; 
    
    [T2dist, T2logbins,nmrName] = loadRawNMRdata(site);
    % load data file
    [d, Dk, T2ML, phi, z, SumEch, kk, lt, lp, SumEch_3s, SumEch_twm, ...
        SumEch_twm_3s] = loadnmrdata2(nmrName); 

    logSumEch = log10(SumEch); 

    %%%%%%%%% Change T2 variable to Sum of Echoes for the inversions. 
    lt = logSumEch; 
    T2ML = SumEch; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Nboot =  1000; % number of bootstrap samples

    
    %[b_boot, n_boot] = bootstrap_fun_mb([lt, kk], Nboot);    % n can vary
    [b_boot, n_boot] = bootstrap_fun([lt, kk], Nboot, n(k));    % n is fixed
   
    meanb(k) = mean(b_boot); 
    meann(k) = n(k); 
    
    medianb(k) = median(b_boot);
    mediann(k) = n(k);

    SOE_K = medianb(k)*(SumEch).^mediann(k);
    k_estimates = [k_estimates SOE_K];
    
    %plotKwithDepth(Dk,z,T2dist,T2logbins,SOE_K,{'DPP','SOE'},{'+'})
    errorEstimate(k) = sum(computeError(Dk, SOE_K));
    
end

plotKestKdpp(Dk,k_estimates,k_names,k_sym)
title(siteName)
 