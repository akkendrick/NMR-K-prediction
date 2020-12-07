% sites_Maurer = {'dpnmr_larned_east','dpnmr_larned_lwph','dpnmr_larned_west',...
%    'dpnmrA11','dpnmrA12','dpnmrC1S','dpnmrC1SE','dpnmrC1SW',...
%    'dpnmr_leque_east','dpnmr_leque_west'};

% % Kansas
%  sites = {'dpnmr_larned_east','dpnmr_larned_lwph','dpnmr_larned_west',...
%     'dpnmrA11','dpnmrA12','dpnmrC1S','dpnmrC1SE','dpnmrC1SW'};

% % Washington
%sites = {'dpnmr_leque_east','dpnmr_leque_west'}

% Wisconsin
sites = {"Site1-WellG5","Site1-WellG6","Site2-WellPN1","Site2-WellPN2"};

SDR_b = [0.0040 0.0024 0.0075 0.0047];
SDR_m = [0 0 0 0];
SDR_n = [2 2 2 2];

SDR_K = @(b,m,n,phi,T2ML) b.*(phi.^m).*(T2ML).^n;

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
    currentSite{kk} = repmat(sites{kk}, length(phi), 1);
    sortedT2dist{kk} = filtT2dist(sortInd,:);
    sortedPhi{kk} = phi(sortInd);
    sortedT2ML{kk} = T2ML(sortInd);
    sortedZ{kk} = z(sortInd);
    sortedLogT2ML{kk} = logT2ML(sortInd);
    sortedLogK{kk} = logK(sortInd);
    sortedSOE{kk} = SumEch(sortInd);
    currentSDRb{kk} = ones(length(phi),1).*SDR_b(kk);
    currentSDRm{kk} = ones(length(phi),1).*SDR_m(kk);
    currentSDRn{kk} = ones(length(phi),1).*SDR_n(kk);
    
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
allSDRb = vertcat(currentSDRb{:});
allSDRm = vertcat(currentSDRm{:});
allSDRn = vertcat(currentSDRn{:});
allSites = vertcat(currentSite{:});


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
allSDRb = allSDRb(allSortInd);
allSDRm = allSDRm(allSortInd);
allSDRn = allSDRn(allSortInd);
allSites = allSites(allSortInd);

allT2linbins = vertcat(siteT2linbins{:});
allT2logbins = vertcat(siteT2logbins{:});

T2Seevers = (allT2ML.^(-1) - T2B^(-1)).^(-1);

T2am = 0;
T2lm = 0;
T2hm = 0;

T2B = 2.2;
thirdT2B = log10(T2B/3);

for jj = 1:length(allZ)
    sliceT2 = allT2dist(jj,:);
    T2linbins = allT2linbins(jj,:);
    T2logbins = allT2logbins(jj,:);
    
    T2logbinsFilt = T2logbins;
    sliceT2Filt = sliceT2;
    
    sliceT2Filt(T2logbins > thirdT2B) = 0;
    T2spacing = diff(T2linbins);
        
%     for j = 2:length(T2linbins)
% 
%         avgSignal = (sliceT2(j)+sliceT2(j-1))/2;
%         signalTemp(j) = avgSignal*T2spacing(j-1);
%         
%     end
%    
%     trapzSummedSignal(jj,1) = trapz(T2linbins', sliceT2');
%     summedSignalCalc(jj,1) = (sum(signalTemp));
%     summedSignalReg(jj,1) = sum(sliceT2);
    
    T2am(jj,1) = sum(T2linbins.*sliceT2)./allPhi(jj); %Arithmetic Mean
    T2lm(jj,1) = 10.^(sum(T2logbins.*sliceT2)./allPhi(jj));
    T2hm(jj,1) = allPhi(jj)./(sum(sliceT2./T2linbins)); 
    
end


%%
figureson = 0;
promRatio = 0;
widthRatio = 0;
largeModeSum = 0;
smallModeSum = 0;

% Compute T2 distribution statistics
for jj = 1:length(allT2ML)

    sliceT2 = allT2dist(jj,:);
    T2linbins = allT2linbins(jj,:);
    T2logbins = allT2logbins(jj,:);
    
    [pks{jj},locs{jj},widths{jj},proms{jj}] = findpeaks(sliceT2);

    if length(pks{jj}) > 1
        T2peakAvg(jj) =  mean(T2linbins(locs{jj}));
        currentProm = proms{jj};
        currentWidth = widths{jj};
        currentLocs = locs{jj};

        promRatio(jj,1) = currentProm(1)/currentProm(2);
        widthRatio(jj,1) = currentWidth(1)/currentWidth(2);

        T2peakMax(jj) = T2linbins(currentLocs(2));
        T2peakMin(jj) = T2linbins(currentLocs(1));
    else
        
        T2peakAvg(jj) = T2linbins(locs{jj});
        promRatio(jj,1) = NaN;
        widthRatio(jj,1) = NaN;
        T2peakMax(jj) = T2linbins(locs{jj});
        T2peakMin(jj) = NaN;
    end

    % Compute summed portions of sliceT2 based on transition between large
    % and small pores
    minLoc = islocalmin(sliceT2);
    [maxVal, maxInd] = max(sliceT2);
    
    stepLocalMax = islocalmax(sliceT2);
    T2localMax = T2linbins(stepLocalMax);
    
    diffPeaks = T2peakMax - T2peakMin;
    
    minInd = find(minLoc == 1);

    % Take minVals less than max
    goodMinLoc = minInd(minInd < maxInd);
    
    % Take the local min with the largest relax time
    goodMinLoc = max(goodMinLoc);
    allGoodMinLoc{jj,1} = goodMinLoc;

    %Fix short/long summed signal! Bin spacing is nonlinear!
    T2spacing = diff(T2linbins);

    if ~isempty(goodMinLoc)
        
        shortSignal(1:goodMinLoc) = sliceT2(1:goodMinLoc); 
        shortSignal(goodMinLoc+1:100) = 0;
        
        [maxVal, indShortMax(jj,1)] = max(shortSignal);
        shortT2Max(jj,1) = T2linbins(indShortMax(jj,1));
        
        longSignal(1:goodMinLoc) = 0; 
        longSignal(goodMinLoc+1:100) = sliceT2(goodMinLoc+1:end);
        
        [maxVal, indLongMax(jj,1)] = max(longSignal);
        longT2Max(jj,1) = T2linbins(indLongMax(jj,1));
        
        longSignals(1:100,jj) = longSignal;
        shortSignals(1:100,jj) = shortSignal;
%         for j = 2:length(T2linbins)
%         
%             avgShort = (shortSignal(j)+shortSignal(j-1))/2; 
%             avgLong = (longSignal(j)+longSignal(j-1))/2; 
% 
%             shortTemp(j) = avgShort*T2spacing(j-1);
%             longTemp(j) = avgLong*T2spacing(j-1);  
%         end
    
        shortSummedSignal(jj) = sum(shortSignal);
        longSummedSignal(jj) = sum(longSignal);
        
        allMinLoc(jj,1) = goodMinLoc;
        allT2min(jj,1) = T2linbins(goodMinLoc);
                
        T2lm_short(jj,1) = 10.^(sum(T2logbins(1:goodMinLoc).*sliceT2(1:goodMinLoc))./shortSummedSignal(jj));
        T2lm_long(jj,1) = 10.^(sum(T2logbins(goodMinLoc+1:end).*sliceT2(goodMinLoc+1:end))./longSummedSignal(jj));
        
        T2_bimodal(jj) = 1;
        
        SDR_K_long(jj) = SDR_K(allSDRb(jj), allSDRm(jj), allSDRn(jj),...
        allPhi(jj), T2lm_long(jj));
    
        SDR_K_short(jj) = SDR_K(allSDRb(jj), allSDRm(jj), allSDRn(jj),...
        allPhi(jj), T2lm_short(jj));

        Sand_frac_approx(jj) = longSummedSignal(jj)/allPhi(jj);
        
        Kparallel(jj) = (shortSummedSignal(jj)/allPhi(jj) * SDR_K_short(jj))+(longSummedSignal(jj)/allPhi(jj) * SDR_K_long(jj));
        Kseries(jj) = 1/(((shortSummedSignal(jj)/allPhi(jj))/SDR_K_short(jj)) + ((longSummedSignal(jj)/allPhi(jj))/SDR_K_long(jj)));
        KregSDR(jj) = SDR_K(allSDRb(jj), allSDRm(jj), allSDRn(jj),...
        allPhi(jj), allT2ML(jj));

        temp = 1/SDR_K_long(jj) + 1/SDR_K_short(jj);
        Kcorr(jj) = 1/temp;
        
        Kreg(jj) = SDR_K(allSDRb(jj), allSDRm(jj), allSDRn(jj),...
        allPhi(jj), allT2ML(jj));
        
        
    else
      
        longSignal(1:100) = sliceT2(1:end);
        longSignals(1:100,jj) = longSignal;
        shortSignals(1:100,jj) = zeros(1,100);
        
%         for j = 2:length(T2linbins)
%             avgLong = (longSignal(j)+longSignal(j-1))/2; 
%             longTemp(j) = avgLong*T2spacing(j-1);  
%         end
        
        shortSummedSignal(jj) = 0;
        longSummedSignal(jj) = sum(longSignal);
        allMinLoc(jj,1) = NaN;
        T2lm_short(jj,1) = NaN;
        allT2min(jj,1) = NaN;
        T2lm_long(jj,1) = T2lm(jj,1);
        
        T2_bimodal(jj) = 0;
        
        [maxVal, indLongMax(jj)] = max(longSignal);
        longT2Max(jj,1) = T2linbins(indLongMax(jj));
        shortT2Max(jj,1) = 0;
        indShortMax(jj,1) = 1;
        
        Sand_frac_approx(jj) = NaN;

        Kparallel(jj) = NaN;
        Kseries(jj) = NaN;
        KregSDR(jj) = NaN;
        
        Kcorr(jj) = SDR_K(allSDRb(jj), allSDRm(jj), allSDRn(jj),...
        allPhi(jj), T2lm_long(jj));
        
        Kreg(jj) = SDR_K(allSDRb(jj), allSDRm(jj), allSDRn(jj),...
        allPhi(jj), allT2ML(jj));

   end

    largeModeSum(jj,1) = longSummedSignal(jj);
    smallModeSum(jj,1) = shortSummedSignal(jj);
    allT2cumsum(jj,:) = cumsum(sliceT2,'reverse');
    
    if indShortMax(jj,1) ~= 1 
        T2PeakDiff(jj,1) = T2linbins(indLongMax(jj,1)) - T2linbins(indShortMax(jj,1));
    else
        T2PeakDiff(jj,1) = 0;
    end    
end

NBoot = 2000;
est = bootstrapK(NBoot, allK, Kparallel);

    
 
  