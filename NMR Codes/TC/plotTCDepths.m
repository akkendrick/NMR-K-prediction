load('optimalCutoffTable_n2_m1_RMSE_2000.mat')

figureson=0;
for kk = 1:length(siteList)
    [T2dist,T2logbins,SEdecayTime,SEdecayUniform,SEdecay,oneDVectors,...
    oneDVectorsUniform, nmrName] = loadAllRawNMRdata(siteList{kk});

    [d, K, T2ML, phi, z, SumEch, log10K, log10T2, log10Porosity,...
    SumEch_3s, SumEch_twm, SumEch_twm_3s] = loadnmrdata2(nmrName);

    phiMatrix{kk} = phi;
    T2MLMatrix{kk} = T2ML;
    KDPPMatrix{kk} = K;
    zMatrix{kk} = z;

    for jj = 1:length(cutoff)
        [K,z,T2dist,T2logbins,kTC_ref_best{jj},bestFitMatrix,totalError,...
        indexQuotient] = computeTCperm_nobtstp(siteList{kk},n,m,median(totalcMatrix(kk,:)),cutoff(jj),figureson);
        
        finalError{kk,jj} = totalError;
        
%         plot(KDPPMatrix{kk},z,'+b','MarkerSize',8)
%         hold on
% 
%         plot(kTC_ref_best{jj},z,'or','MarkerSize',8)
%         grid on
%         box on
%         
%         set(gca,'XScale','log')
%         fileString = strcat('cutoff= ', string(cutoff(jj)),' ',siteList{kk},'.png');
%         title(fileString)
%         set(gca, 'YDir','reverse')
%         print('-dpng','-r300',fileString)
%         close
    end
end
  