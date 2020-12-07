%% Script to estimate b at each site using grid search

clear

%all_names = {'A1', 'C1', 'dpnmr_larned_east', 'dpnmr_larned_west', ...
%  'dpnmr_larned_lwph', 'dpnmr_leque_east', 'dpnmr_leque_west'}; 

all_names = {'Site1-WellG6','Site1-WellG5', 'Site2-WellPN1', 'Site2-WellPN2'};


baseDir = '/Volumes/GoogleDrive/My Drive/Stanford/USGS Project/Field Data/USGS Data/';


figureson = 1;

for k = 1:length(all_names)
    name = all_names{k}
    site = name;
    if strcmp(site,'Site1-WellG5')
            name = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1';
            nmrName = name;

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
    
    
    %name = 'all_data'
    [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(name); 
    
    logSumEch = log10(SumEch); 
    
    if figureson == 1
        n = 2;      % fix n for the grid search, allow m to vary
        [bs, ms, r] = grid_search(logT2ML, logPhi, logK, n);
    end

    bs_grid(:,k) = bs;
    ms_grid(:,k) = ms;
    

end

%%
% Print results for m = 4 
b = 10.^(bs_grid(4,:))