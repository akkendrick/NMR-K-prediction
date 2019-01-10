%% Script to estimate b at each site using grid search

clear

all_names = {'A1', 'C1', 'dpnmr_larned_east', 'dpnmr_larned_west', ...
  'dpnmr_larned_lwph', 'dpnmr_leque_east', 'dpnmr_leque_west'}; 

figureson = 1;

for k = 1:length(all_names)
    name = all_names{k}
    
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