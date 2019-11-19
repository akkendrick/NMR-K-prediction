% Range over pairs of m and n values
close all
clear

load enso 

figureson = 0;


% Kansas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  siteList = [{'dpnmr_larned_east'} {'dpnmr_larned_lwph'} {'dpnmr_larned_west'}... 
%      {'dpnmrA11'} {'dpnmrA12'} {'dpnmrC1S'} {'dpnmrC1SE'} {'dpnmrC1SW'}];
% 
% refCutoff = [33 33 33 33 33 33 33 33];
% idealCutoff = [84 20 46 40 52 164 174 92];
% 
% SDR_b = [0.027 0.016 0.024 0.032 0.032 0.028 0.028 0.028];
% SDR_m = [0 0 0 0 0 0 0 0];
% SDR_n = [2 2 2 2 2 2 2 2];
  
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Washington
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
siteList = [{'dpnmr_leque_east'}];% {'dpnmr_leque_west'}];

refCutoff = [33 33]./1000;
idealCutoff = [98 106]./1000;

SDR_b = [0.036 0.022];
SDR_m = [0 0];
SDR_n = [2 2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wisc
% siteList = [{'Site1-WellG5'} {'Site1-WellG6'},...
%     {'Site2-WellPN1'} {'Site2-WellPN2'} ];
% 
% refCutoff = [33, 33, 33, 33]./1000;
% idealCutoff = [1060, 548, 226, 318]./1000;
% 
% %idealCutoff = [1000 1000 1000 1000]./1000;
% 
% SDR_b = [0.0024 0.0040 0.0075 0.0047];
% SDR_m = [0 0 0 0];
% SDR_n = [2 2 2 2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SDR_K = @(b,m,n,phi,T2ML) b.*(phi.^m).*(T2ML).^n;

for kk = 1:length(siteList)
   
    [T2dist,T2logbins,SEdecayTime,SEdecayUniform,SEdecay,oneDVectors,...
    oneDVectorsUniform, nmrName] = loadAllRawNMRdata(siteList{kk});

    [d, K, T2ML, phi, z, SumEch, log10K, log10T2, log10Porosity,...
        SumEch_3s, SumEch_twm, SumEch_twm_3s] = loadnmrdata2(nmrName);
    
    phiMatrix{kk} = phi;
    T2MLMatrix{kk} = T2ML;
    KMatrix{kk} = K;
    zMatrix{kk} = z;
    
    
    [K,z,T2dist,T2logbins,kTC_ref_best{kk},bestFitMatrix,totalError,...
        indexQuotient] = computeTCperm_2(siteList{kk},SDR_n(kk),SDR_m(kk),refCutoff(kk),figureson);
    
    [K,z,T2dist,T2logbins,kTC_ideal_best{kk},bestFitMatrix,totalError,...
    indexQuotient] = computeTCperm_2(siteList{kk},SDR_n(kk),SDR_m(kk),idealCutoff(kk),figureson);

    SDR_Kest{kk} = SDR_K(SDR_b(kk),SDR_m(kk),SDR_n(kk),phiMatrix{kk},T2MLMatrix{kk}); 
    
    SDRError{kk} = estimateKdiffFactor(KMatrix{kk},SDR_Kest{kk},1);  
    TC_refError{kk} = estimateKdiffFactor(KMatrix{kk},kTC_ref_best{kk},1);  
    TC_idealError{kk} = estimateKdiffFactor(KMatrix{kk},kTC_ideal_best{kk},1);  
    
    meanSDRError(kk) = mean(SDRError{kk});
    meanTC_refError(kk) = mean(TC_refError{kk});
    meanTC_idealError(kk) = mean(TC_idealError{kk});
       
    
end