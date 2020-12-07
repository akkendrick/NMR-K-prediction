% Find ideal T2 with constant b 

% % Wisc optimized
% SDR_b = 0.0043;

% Maurer and Knight
SDR_b = 0.036;

% % Wisconsin
sites = {'Site1-WellG5','Site1-WellG6','Site2-WellPN1','Site2-WellPN2'};

SDR_m = 0;
SDR_n = 2;

SDR_K = @(b,m,n,phi,T2ML) b.*(phi.^m).*(T2ML).^n;
SDR_T2 = @(b,m,n,phi,K) (K./(b.*(phi.^m))).^(1/n);


for kk = 1:length(sites)   
    
    [T2dist, T2logbins, nmrName] = loadRawNMRdata(sites{kk});
    
   % gammaEMIdepth = [];
   % [gammaEMIdepth, gamma, EMI] = loadGammaEMIData(sites{kk});
    
    [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(nmrName);

    T2depths = T2dist(:,1);
    T2data = T2dist(:,2:end);
    
    T2axis = 10.^T2logbins;
    T2axis = T2axis * 10^3;
        
    T2B = 2.0045;
    
    for jj = 1:length(z)
       [minVal, index] = min(abs(T2depths - z(jj)));
       filtIndices(jj) = index; 
       filtT2dist(jj,:) = T2data(index,:);        
       filtT2depths(jj) = T2depths(index);
    end
    
    [sortedK{kk}, sortInd] = sort(K);
    sortedT2dist{kk} = filtT2dist(sortInd,:);
    sortedPhi{kk} = phi(sortInd);
    sortedT2ML{kk} = T2ML(sortInd);
    sortedZ{kk} = z(sortInd);
    sortedLogT2ML{kk} = logT2ML(sortInd);
    sortedLogK{kk} = logK(sortInd);
    sortedSOE{kk} = SumEch(sortInd);
    
    T2logbins = ones(length(phi),length(T2logbins)).*T2logbins;
    
    siteT2linbins{kk} = 10.^(T2logbins);
    siteT2logbins{kk} = T2logbins;
end

allK = vertcat(sortedK{:});
allT2dist = vertcat(sortedT2dist{:});
allPhi = vertcat(sortedPhi{:});
allT2ML = vertcat(sortedT2ML{:});
allZ = vertcat(sortedZ{:});
allLogT2ML = vertcat(sortedLogT2ML{:});
allLogK = vertcat(sortedLogK{:});
allSOE = vertcat(sortedSOE{:});


allT2linbins = vertcat(siteT2linbins{:});
allT2logbins = vertcat(siteT2logbins{:});

%[allK, allSortInd] = sort(allK);
[sortedAllK, allSortInd] = sort(allK);

allK = allK(allSortInd);
allT2dist = allT2dist(allSortInd,:);
allPhi = allPhi(allSortInd);
allT2ML = allT2ML(allSortInd);
allZ = allZ(allSortInd);
allLogT2ML = allLogT2ML(allSortInd);
allLogK = allLogK(allSortInd);
allSOE = allSOE(allSortInd);

figureson = 0;

for jj = 1:length(allZ)
    sliceT2 = allT2dist(jj,:);
    T2linbins = allT2linbins(jj,:);
    T2logbins = allT2logbins(jj,:);

    T2_ideal(jj) = SDR_T2(SDR_b,SDR_m,SDR_n,allPhi(jj),allK(jj));
    K_SDR_ideal(jj) = SDR_K(SDR_b,SDR_m,SDR_n,allPhi(jj),T2_ideal(jj));
    
    if figureson == 1

        [val, index] = min(abs(allT2ML(jj) - T2linbins));
        [val, index2] = min(abs(T2_ideal(jj) - T2linbins));

        plot(T2linbins, sliceT2, 'LineWidth',2)

        grid on
        box on
        hold on

        plot(T2linbins(index), sliceT2(index), 'k*','MarkerSize',10)
        plot(T2linbins(index2), sliceT2(index2), 'r*','MarkerSize',10)

        hold off
        set(gca,'XScale','log')

        legendString = strcat('K= ', sprintf('%d',allK(jj)),'m/s');
        fileString = strcat('K= ', sprintf('%d',allK(jj)),'.png');

        legend(legendString,'T2ML','Ideal SDR T2ML','Location','northwest')
        print('-dpng','-r300',fileString)
    end
    
end


