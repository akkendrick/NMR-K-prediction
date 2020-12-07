% all_names = {'A1','C1', 'dpnmr_leque_east', 'dpnmr_leque_west', 'dpnmr_larned_east', ...
%    'dpnmr_larned_west', 'dpnmr_larned_lwph', 'gems_all', 'dpnmr_leque_all','dpnmr_larned_all',  'all_data'}; 

%all_names = {'Site1-WellG5','Site1-WellG5above','Site1-WellG5below','Site1-WellG6','Site1-WellG6above','Site1-WellG6below','Site2-WellPN1','Site2-WellPN2'};
all_names = [{'Site1-WellG6'} {'Site1-WellG5'} {'Site2-WellPN1'} {'Site2-WellPN2'}];

baseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/USGS Data/';
%baseDir = 'I:\My Drive\USGS Project\USGS Data\';


C = {'-*k', '-*r', '-.b', '-.g', '-oc', '-om', '-or', '-*b', '-.b', '-ob', '--k'}; 
figure; 
hold on

for k = 1:length(all_names)
     clearvars -except all_names k rnorm C baseDir
     site = all_names{k}
     
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

    %% Test which SOE works best
    % 1 - SDR
    % 2 - SOE
    % 3 - SOE - 3 seconds (extrapolated)
    % 4 - SOE - time-weighted mean
    % 5 - SOE - 3 second time-weighted mean

    vars = [lt, log10(SumEch)];%, log10(SumEch_3s), log10(SumEch_twm), log10(SumEch_twm_3s)]; 
    numvars = size(vars, 2);
    rnorm = []; 
    for i = 1:numvars
        var = vars(:,i); 
        r = kk - (var\kk)*var; 
        rnorm(i) = norm(r)./length(r); 
    end
    
    plot([1:numvars], rnorm, C{k})

end
set(gca, 'xtick', [1:5])
labels = {'T_{2ML}', 'SOE', 'SOE - 3s', 'SOE - TWM', 'SOE - 3s, TWM'}; 
set(gca, 'xticklabel', labels)
xlim([0.5, 5.5])
ax = gca;
ax.XTickLabelRotation = 45;

%names = {'GEMS2 - A1', 'GEMS2 - C1', 'Leque East', 'Leque West', 'Larned East', 'Larned West', 'Larned C', 'GEMS2 - all', 'Leque - all', 'Larned - all', 'All Data'}; 
names = {'G5','G5a','G5b','G6','G6a','G6b','PN1','PN2'}
legend(names')


%% Now look at SOE
Nboot =  2000; % number of bootstrap samples

for k = 1:length(all_names)

     site = all_names{k}
    
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
        name = site;
        nmrName = site;
    end



    % load data file
    [d, Dk, T2ML, phi, z, SumEch, kk, lt, lp, SumEch_3s, SumEch_twm, ...
        SumEch_twm_3s] = loadnmrdata2(name); 
    logSumEch = log10(SumEch); 
%     logSumEch_3s = log10(SumEch_3s); 
%     logSumEch_twm = log10(SumEch_twm); 
%     logSumEch_twm_3s = log10(SumEch_twm_3s); 
    %%%%%%%%% Change T2 variable to Sum of Echoes for the inversions. 
    lt = logSumEch; 
    T2ML = SumEch; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    n = 2;
    %[b_boot, n_boot] = bootstrap_fun_mb([lt, kk], Nboot);    % n can vary
    [b_boot, n_boot] = bootstrap_fun([lt, kk], Nboot, n);    % fix n
    meanb(k) = mean(b_boot); 
    meann(k) = mean(n_boot); 
    
    medianb(k) = median(b_boot);
    mediann(k) = median(n_boot);

    SOE_K = medianb(k)*(SumEch).^mediann(k);
    
    plotKwithDepth(Dk,z,T2dist,T2logbins,SOE_K,{'DPP','SOE'},{'+'})
    errorEstimate(k) = sum(computeError(Dk, SOE_K));
    
   % graph_correlations([b_boot, n_boot], 2, {'log_{10}(b)', 'n'}, 1, 0)
end

figure; 
plot([1:k], meann)
