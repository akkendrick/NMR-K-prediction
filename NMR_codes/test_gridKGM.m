clear

name = 'A1'
%name = 'A1'
%name = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1';

[d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(name); 

[bestTau, bestRho, r] = grid_search_kgm(logT2ML, logPhi, K);

[bestTau, bestRho, r] = grid_search_KGM_akk(T2ML, phi, K);