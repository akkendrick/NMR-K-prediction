load('SDR_wiscAll_noFilt_m0_n2.mat')

indSiteList = {'Site1-WellG5','Site1-WellG6','Site2-WellPN1','Site2-WellPN2'};

for kk=1:length(indSiteList)

    site = indSiteList{kk}
    
    if strcmp(site,'Site1-WellG5')
        name = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1'; 
        nmrName = name;
        saveName = 'G5';

    elseif strcmp(site,'Site1-WellG6')
        name = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';
        nmrName = name;
        saveName = 'G6';

    elseif strcmp(site,'Site2-WellPN1')
        name = 'Pl_W1_Tr5_20x_MPp75aLS_F1n2_wRIN_wRFI_Reg50_Va1';
        nmrName = name;
        saveName = 'PN1';

    elseif strcmp(site,'Site2-WellPN2')
        name = 'W2_Tr5_20x_MPp75aLS_Reg50_wRIN_wRFI_Va1';
        nmrName = name;
        saveName =                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      'PN2';

    else
        name = site;
        nmrName = site;
        saveName = 'wisc_all';
    end
    
    mVal = m;
    nVal = n;
    bVal = totalbMatrix(1);
    
    [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
        SumEch_twm_3s] = loadnmrdata2(nmrName); 
    
     SDR_K = @(b,m,n,phi,T2ML) b.*(phi.^m).*(T2ML).^n;

     SDR_Kest{:,kk} = SDR_K(bVal, mVal, nVal, phi, T2ML);
     SDRErrorFactor{:,kk} = mean(estimateKdiffFactor(K,SDR_Kest{:,kk},1));
     
    
     
end

