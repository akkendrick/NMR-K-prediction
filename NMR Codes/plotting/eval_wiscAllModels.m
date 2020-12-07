clear
close all

load 'SDR_wiscAll_noFilt_mRange_nRange.mat'

site = 'wisc_all';

[d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(site); 

 SDR_K = @(b,m,n,phi,T2ML) b.*(phi.^m).*(T2ML).^n;

 m = [0 1 2 4 8 0 1 2 4 8];
 n = [2 2 2 2 2 1 1 1 1 1];
 
 for kk = 1:length(totalbMatrix)
     
     mVal = m(kk);
     nVal = n(kk);
     bVal = totalbMatrix(1,kk);
     
     SDR_Kest{:,kk} = SDR_K(bVal, mVal, nVal, phi, T2ML);
     
     SDRErrorFactor{:,kk} = mean(estimateKdiffFactor(K,SDR_Kest{:,kk},1));
     
 end
 
 plotKestKdpp_wdepth(K,[SDR_Kest{1} SDR_Kest{9}],[SDR_Kest{6} SDR_Kest{9}],z)
 