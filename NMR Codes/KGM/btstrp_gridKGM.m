% Bootstrap KGM model

% siteList = [{'Site1-WellG5'} {'Site1-WellG6'}  {'Site2-WellPN1'} {'Site2-WellPN2'}];

siteList = {'dpnmr_larned_east','dpnmr_larned_lwph','dpnmr_larned_west',...
   'dpnmrA11','dpnmrA12','dpnmrC1S','dpnmrC1SE','dpnmrC1SW',...
   'dpnmr_leque_east','dpnmr_leque_west'};

nBoot = 2000;

% Make sure to set temp to 10.6 degrees to get close to Wisc data

parfor kk = 1:length(siteList)
   [T2dist, T2logbins, nmrName] = loadRawNMRdata(siteList{kk});
   
    [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(nmrName); 
    
    m = 1;

    data{kk} = [logT2ML phi logK];
       
    bootKGM{kk} = bootstrp(nBoot,@(X) btstrp_kgm(X,m), data{kk});
    
%    [KGM_lk, bestTau, bestRho, r] = grid_search_kgm(logT2ML, phi, logK, m);

    %test = btstrp_kgm(data, m);


%     KGM_k = 10.^KGM_lk;
%     display('RMSE Error')
%     bestError = computeError(K, KGM_k)

    
end