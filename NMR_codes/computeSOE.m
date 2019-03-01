% computeSOE

% Analyze SOE by plotting K_SOE vs K_DPP
baseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/USGS Data/';

n = [2 1 0];

k_estimates = [];
k_names = {'1:1','SOE n=2','SOE n=1','SOE n=0'};
k_sym = {'+','*','o'};

sites = {'


for k = 1:length(n)
    k
    site = 'Site2-WellPN2'
    siteName = site; 
    
    if strcmp(site,'Site1-WellG5')
            name = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1';
            nmrName = name;

            in1 = [baseDir site '/' name '/' name '_T2_dist' '.txt']; 
            in2 = [baseDir site '/' name '/' name '_T2_bins_log10s' '.txt']; 

            T2dist = load(in1); 
            T2logbins = load(in2);
        elseif  strcmp(site,'Site1-WellG5above')
            site = 'Site1-WellG5';
            name = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1';
            nmrName = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1_above';

            in1 = [baseDir site '/' name '/' name '_T2_dist' '.txt']; 
            in2 = [baseDir site '/' name '/' name '_T2_bins_log10s' '.txt']; 

            T2dist = load(in1); 
            T2logbins = load(in2);
        elseif  strcmp(site,'Site1-WellG5below')
            site = 'Site1-WellG5';
            name = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1';
            nmrName = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1_below';

            in1 = [baseDir site '/' name '/' name '_T2_dist' '.txt']; 
            in2 = [baseDir site '/' name '/' name '_T2_bins_log10s' '.txt']; 

            T2dist = load(in1); 
            T2logbins = load(in2);
        elseif strcmp(site,'Site1-WellG6')
            name = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';
            nmrName = name;

            in1 = [baseDir site '/' name '/' name '_T2_dist' '.txt']; 
            in2 = [baseDir site '/' name '/' name '_T2_bins_log10s' '.txt']; 

            T2dist = load(in1); 
            T2logbins = load(in2);
        elseif strcmp(site,'Site1-WellG6above')
            site = 'Site1-WellG6';
            name = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';
            nmrName = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1_above';

            in1 = [baseDir site '/' name '/' name '_T2_dist' '.txt']; 
            in2 = [baseDir site '/' name '/' name '_T2_bins_log10s' '.txt']; 

            T2dist = load(in1); 
            T2logbins = load(in2);
        elseif strcmp(site,'Site1-WellG6below')
            site = 'Site1-WellG6';
            name = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';
            nmrName = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1_below';
            in1 = [baseDir site '/' name '/' name '_T2_dist' '.txt']; 
            in2 = [baseDir site '/' name '/' name '_T2_bins_log10s' '.txt']; 

            T2dist = load(in1); 
            T2logbins = load(in2);
        elseif strcmp(site,'Site2-WellPN1')
            name = 'Pl_W1_Tr5_20x_MPp75aLS_F1n2_wRIN_wRFI_Reg50_Va1';
            nmrName = name;
            in1 = [baseDir site '/' name '/' name '_T2_dist' '.txt']; 
            in2 = [baseDir site '/' name '/' name '_T2_bins_log10s' '.txt']; 

            T2dist = load(in1); 
            T2logbins = load(in2);
        elseif strcmp(site,'Site2-WellPN2')
            name = 'W2_Tr5_20x_MPp75aLS_Reg50_wRIN_wRFI_Va1';
            nmrName = name;

            in1 = [baseDir site '/' name '/' name '_T2_dist' '.txt']; 
            in2 = [baseDir site '/' name '/' name '_T2_bins_log10s' '.txt']; 

            T2dist = load(in1); 
            T2logbins = load(in2);
        else
            nmrName = site;
    end

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
    meann(k) = mean(n_boot); 
    
    medianb(k) = median(b_boot);
    mediann(k) = median(n_boot);

    SOE_K = medianb(k)*(SumEch).^mediann(k);
    k_estimates = [k_estimates SOE_K];
    
    %plotKwithDepth(Dk,z,T2dist,T2logbins,SOE_K,{'DPP','SOE'},{'+'})
    errorEstimate(k) = sum(computeError(Dk, SOE_K));
    
end

plotKestKdpp(Dk,k_estimates,k_names,k_sym)
title(siteName)
