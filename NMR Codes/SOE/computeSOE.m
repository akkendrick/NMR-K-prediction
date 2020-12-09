% computeSOE

% Analyze SOE by plotting K_SOE vs K_DPP
%baseDir = '/Volumes/GoogleDrive/My Drive/Stanford/USGS Project/Field Data/USGS Data/';
baseDir = 'I:\My Drive\Stanford\USGS Project\Field Data\USGS Data\';

n = 1;

k_estimates = [];
k_names = {'1:1','SOE n=2','SOE n=1','SOE n=0'};
k_sym = {'+','*','o'};

siteList = [{'Site1-WellG5'} {'Site1-WellG6'}  {'Site2-WellPN1'} {'Site2-WellPN2'}];

% siteList = {'dpnmr_larned_east','dpnmr_larned_lwph','dpnmr_larned_west',...
%   'dpnmrA11','dpnmrA12','dpnmrC1S','dpnmrC1SE','dpnmrC1SW',...
%   'dpnmr_leque_east','dpnmr_leque_west'};


for k = 1:length(siteList)

    site = siteList{k}; 
    siteName = site
    
    [T2dist, T2logbins,nmrName] = loadRawNMRdata(site);

    % load data file
    [d, Dk, T2ML, phi, z, SumEch, kk, lt, lp, SumEch_3s, SumEch_twm, ...
        SumEch_twm_3s] = loadnmrdata2(nmrName); 

    logSumEch = log10(SumEch); 

    %%%%%%%%% Change T2 variable to Sum of Echoes for the inversions. 
    lt = logSumEch; 
    T2ML = SumEch; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Nboot =  2000; % number of bootstrap samples

    
    %[b_boot, n_boot] = bootstrap_fun_mb([lt, kk], Nboot);    % n can vary
    [b_boot, n_boot] = bootstrap_fun([lt, kk], Nboot, n);    % n is fixed
   
    meanb(k) = mean(b_boot); 
    meann(k) = mean(n_boot); 
    
    medianb(k) = median(b_boot)
    mediann(k) = median(n_boot)

    b_boot_all{k} = b_boot;
    
    SOE_K = medianb(k)*(SumEch).^mediann(k);
    k_estimates{k} = SOE_K;
    
    %plotKwithDepth(Dk,z,T2dist,T2logbins,SOE_K,{'DPP','SOE'},{'+'})
    errorEstimate(k) = sum(computeError(Dk, SOE_K));
    
end

save('SOE_dat.mat','siteList','b_boot_all','medianb','mediann','k_estimates','errorEstimate')
% plotKestKdpp(Dk,k_estimates,k_names,k_sym)
% title(siteName)
