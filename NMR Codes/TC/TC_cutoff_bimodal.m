% Look for avg cutoff from bimodal T2 distributions

sites = {'Site1-WellG5','Site1-WellG6','Site2-WellPN1','Site2-WellPN2'};


for kk = 1:length(sites)   
    
    site = squeeze(sites{kk})
    
    [T2dist, T2logbins, nmrName] = loadRawNMRdata(sites{kk});
    
    [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(site);

    T2depths = T2dist(:,1);
    T2data = T2dist(:,2:end);
        
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


[sortedAllK, allSortInd] = sort(allK);

allK = allK(allSortInd);
allT2dist = allT2dist(allSortInd,:);
allPhi = allPhi(allSortInd);
allT2ML = allT2ML(allSortInd);
allZ = allZ(allSortInd);
allLogT2ML = allLogT2ML(allSortInd);
allLogK = allLogK(allSortInd);
allSOE = allSOE(allSortInd);

allT2linbins = vertcat(siteT2linbins{:});
allT2logbins = vertcat(siteT2logbins{:});

for jj = 1:length(allZ)
    sliceT2 = allT2dist(jj,:);
    T2linbins = allT2linbins(jj,:);
    T2logbins = allT2logbins(jj,:);
    
    T2logbinsFilt = T2logbins;
    sliceT2Filt = sliceT2;
    
    T2spacing = diff(T2linbins);
    
    T2am(jj,1) = sum(T2linbins.*sliceT2)./allPhi(jj); %Arithmetic Mean
    T2lm(jj,1) = 10.^(sum(T2logbins.*sliceT2)./allPhi(jj));
    
    %T2lm_test(jj,1) = (sum(T2linbins.*sliceT2)./allPhi(jj));

    T2hm(jj,1) = allPhi(jj)./(sum(sliceT2./T2linbins)); 
    
end

figureson = 0;
promRatio = 0;
widthRatio = 0;
largeModeSum = 0;
smallModeSum = 0;

% Compute T2 distribution statistics
for jj = 1:length(allT2ML)
%for jj = 2:2

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
%         for j = 2:length(T2linbins)
%         
%             avgShort = (shortSignal(j)+shortSignal(j-1))/2; 
%             avgLong = (longSignal(j)+longSignal(j-1))/2; 
% 
%             shortTemp(j) = avgShort*T2spacing(j-1);
%             longTemp(j) = avgLong*T2spacing(j-1);  
%         end
    

        is_bimodal(jj) = 1;
        bimodalT2(jj,1) = T2linbins(goodMinLoc);

    else
              
        is_bimodal(jj) = 0;
        bimodalT2(jj,1) = NaN;

    end
   
end

disp('Median location of min between T2 modes in s')
median(bimodalT2,'omitnan')



