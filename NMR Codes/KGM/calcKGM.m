

sites = ["Site1-WellG5","Site1-WellG6","Site2-WellPN1","Site2-WellPN2"];

m = 1;

for k=1:length(sites)
    site = sites{k}
    
    [T2dist, T2logbins,nmrName] = loadRawNMRdata(site);
    
    % load data file
    [d, Dk, T2ML, phi, z, SumEch, kk, lt, lp, SumEch_3s, SumEch_twm, ...
        SumEch_twm_3s] = loadnmrdata2(nmrName);
    
    DPP_K{k} = Dk;
    logT2ML = log10(T2ML);
    logK = log10(Dk);
    
    [KGM_lk{k}, bestTau{k}, bestRho{k}, r{k}] = grid_search_kgm(logT2ML, phi, logK, m);
    
    KGM_K{k} = 10.^KGM_lk{k};
    
    bestError{k} = computeError(DPP_K{k}, KGM_K{k});
    [errorSign{k}, errorFactor{k}] = estimateKdiffFactor_withSign(DPP_K{k},KGM_K{k},1);

end
