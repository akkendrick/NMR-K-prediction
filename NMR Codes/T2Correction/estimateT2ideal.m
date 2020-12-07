clear


% Kansas
%sites = {'dpnmr_larned_east','dpnmr_larned_lwph','dpnmr_larned_west',...
%   'dpnmrA11','dpnmrA12','dpnmrC1S','dpnmrC1SE','dpnmrC1SW'};



% Washington
%sites = {'dpnmr_leque_east','dpnmr_leque_west'}
% 

% Maurer and Knight
sites = {'dpnmr_larned_east','dpnmr_larned_lwph','dpnmr_larned_west',...
  'dpnmrA11','dpnmrA12','dpnmrC1S','dpnmrC1SE','dpnmrC1SW',...
  'dpnmr_leque_east','dpnmr_leque_west'};



% Wisconsin
% sites = {'Site1-WellG5','Site1-WellG6','Site2-WellPN1','Site2-WellPN2'};

 
 % all sites
%  sites = {'Site1-WellG5','Site1-WellG6','Site2-WellPN1','Site2-WellPN2',...
%    'dpnmr_larned_east','dpnmr_larned_lwph','dpnmr_larned_west',...
%    'dpnmr_leque_east','dpnmr_leque_west','dpnmrA11','dpnmrA12',...
%    'dpnmrC1S','dpnmrC1SE','dpnmrC1SW'};

for kk = 1:length(sites)   
    
    [T2dist, T2logbins, nmrName] = loadRawNMRdata(sites{kk});
    [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(nmrName);

    T2depths = T2dist(:,1);
    T2data = T2dist(:,2:end);
    
    T2axis = 10.^T2logbins;
    T2axis = T2axis * 10^3;
        
    T2B = 2.2;
    
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
    
     T2logbins = ones(length(phi),length(T2logbins)).*T2logbins;
    
    siteT2linbins{kk} = 10.^(T2logbins);
    siteT2logbins{kk} = T2logbins;
    
    currentLogT2ML = sortedLogT2ML{kk};
    currentK = sortedK{kk};
    currentLogK = sortedLogK{kk};
    
    T2MLlogCutoff = -0.78;
    smallLogK = currentLogK(currentLogT2ML < T2MLlogCutoff);
    largeLogK = currentLogK(currentLogT2ML > T2MLlogCutoff);
    smallLogT2ML = currentLogT2ML(currentLogT2ML < T2MLlogCutoff);
    largeLogT2ML = currentLogT2ML(currentLogT2ML > T2MLlogCutoff);

    %corrFactor = 0.4;
    corrFactor = 0;
    
    % CORR LOG T2ML HAS NOT BEEN SORTED PROPERLY
    
    corrLargeLogT2ML = largeLogT2ML - corrFactor;
    corrlogT2ML{kk} = [smallLogT2ML; corrLargeLogT2ML];
    corrlogK{kk} = [smallLogK; largeLogK];
end

allK = vertcat(sortedK{:});
allT2dist = vertcat(sortedT2dist{:});
allPhi = vertcat(sortedPhi{:});
allT2ML = vertcat(sortedT2ML{:});
allZ = vertcat(sortedZ{:});
allLogT2ML = vertcat(sortedLogT2ML{:});
allLogK = vertcat(sortedLogK{:});
allCorrLogK = vertcat(corrlogK{:});
%allGamma = vertcat(sortedGamma{:});
%allGammaDepths = vertcat(sortedGammaDepths{:});
%allEMI = vertcat(sortedEMI{:});

allT2linbins = vertcat(siteT2linbins{:});
allT2logbins = vertcat(siteT2logbins{:});
allCorrT2ML = 10.^(vertcat(corrlogT2ML{:}));
allCorrK = 10.^(vertcat(corrlogK{:}));
T2Seevers = (allT2ML.^(-1) - T2B^(-1)).^(-1);

%%
T2am = 0;
T2lm = 0;
T2hm = 0;

T2B = 1.5;
thirdT2B = log10(T2B/3);

for jj = 1:length(allZ)
    sliceT2 = allT2dist(jj,:);
    T2linbins = allT2linbins(jj,:);
    T2logbins = allT2logbins(jj,:);
    
    T2logbinsFilt = T2logbins;
    sliceT2Filt = sliceT2;
    
    sliceT2Filt(T2logbins > thirdT2B) = 0;

    
    T2am(jj,1) = sum(T2linbins.*sliceT2)./allPhi(jj); %Arithmetic Mean
    T2lm(jj,1) = 10.^(sum(T2logbins.*sliceT2)./allPhi(jj));
    %T2lm_filt(jj,1) = 10.^(sum(T2logbinsFilt.*sliceT2Filt)./allPhi(jj));
    T2hm(jj,1) = allPhi(jj)./(sum(sliceT2./T2linbins)); 
    
    % Attempt to discritize T2 data
    
    [T2discStep, T2discEdgeStep] = discretize(sliceT2, 20);
    T2disc(jj,:) = T2discStep;
    
    discSteps = 15;
    T2binNum = discretize(T2logbins, discSteps);
    for kk = 1:discSteps
        sliceT2disc(kk) = sum(sliceT2(T2binNum == kk));
    end
    
    allSliceT2disc(jj,:) = sliceT2disc;
    
end

%%
promRatio = 0;
widthRatio = 0;
largeModeSum = 0;
smallModeSum = 0;
% Compute T2 distribution statistics
for jj = 1:length(allT2ML)
    sliceT2 = allT2dist(jj,:);
    T2linbins = allT2linbins(jj,:);
    T2logbins = allT2logbins(jj,:);

    
    % Try removing long time data to see what happens
    longTimeCutoff = 0.15; %s
    [minVal, cutoffBin] = min(abs(T2linbins - longTimeCutoff));
    
    filteredT2Dist = allT2dist(:,1:cutoffBin);
    filteredSlice = filteredT2Dist(jj,:);
    filteredT2ml(jj,1) = 10.^(sum(T2logbins(1:cutoffBin).*filteredSlice)./sum(filteredSlice));
    
    
    [pks{jj}, locs{jj},widths{jj},proms{jj}] = findpeaks(sliceT2);

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
    maxRelaxTime(jj, 1) = T2linbins(maxInd);
    minInd = find(minLoc == 1);

    % Take minVals less than max
    goodMinLoc = minInd(minInd < maxInd);
    
    % Take the local min with the largest relax time
    goodMinLoc = max(goodMinLoc);
    
    goodMinLoc = cutoffBin;
    

    if ~isempty(goodMinLoc)
        shortSummedSignal(jj) = sum(sliceT2(1:goodMinLoc));
        longSummedSignal(jj) = sum(sliceT2(goodMinLoc:end));
        allMinLoc(jj,1) = goodMinLoc;
        minModeT2(jj,1) = T2linbins(allMinLoc(jj,1));

        
        T2lm_short(jj,1) = 10.^(sum(T2logbins(1:goodMinLoc).*sliceT2(1:goodMinLoc))./shortSummedSignal(jj));
        T2lm_long(jj,1) = 10.^(sum(T2logbins(goodMinLoc:end).*sliceT2(goodMinLoc:end))./longSummedSignal(jj));

    else
        shortSummedSignal(jj) = 0;
        longSummedSignal(jj) = sum(sliceT2);
        allMinLoc(jj,1) = NaN;
        T2lm_short(jj,1) = NaN;
        T2lm_long(jj,1) = NaN;
        minModeT2(jj,1) = NaN;


    end

    largeModeSum(jj,1) = longSummedSignal(jj);
    smallModeSum(jj,1) = shortSummedSignal(jj);
    
    allT2cumsum(jj,:) = cumsum(sliceT2,'reverse');

    maxModeT2(jj,1) = T2linbins(maxInd);
end

%%
%optimalbAll = 0.0043;
%optimalbEqualWeight = 0.0020;
%optimalbAll = 3*10^-2;

% Wisconsin Optimal b
%optimalb = 0.000773852; %Lowest K b
%optimalb = 0.0048 % mid range K b
%optimalb = 0.006 % high range K b

% Maurer and Knight Optimal b
%optimalb = 0.0011 % Lowest K b
%optimalb = 0.0259 % mid range K b 
%optimalb = 0.0596 % high range K b
optimalb = 0.0278 % Opt value for all sites in paper 

SDR_b = @(K,m,n,phi,T2ML) K./((phi.^m).*(T2ML).^n);
SDR_T2 = @(K,m,n,phi,b) (K./((phi.^m).*(b))).^(1/n);
SDR_K = @(b,m,n,phi,T2ML) b.*(phi.^m).*(T2ML).^n;

m = 0;
n = 2;

%allbProfile = SDR_b(allK,m,n,allPhi,allT2ML);
allbProfile = SDR_b(allK,m,n,allPhi,allT2ML);
allT2Profile = SDR_T2(allK,m,n,allPhi,optimalb);
allSDRKProfile = SDR_K(optimalb,m,n,allPhi,allT2ML);
%%

T2diff = allT2ML - allT2Profile;
T2diff(T2diff < 0) = NaN;

figure(4)
hold on
grid on
box on
scatter(allK, T2diff, 30, 'Filled')
scatter(allK, T2lm_short, 30, 'Filled')
set(gca,'XScale','log')
set(gca,'YScale','log')

xlabel('K (m/s)')
ylabel('T2 Estimator (s)')

%%
figure(5)
grid on
box on
hold on
scatter(T2diff, maxModeT2,30,'Filled')

%%

% figure(2)
% grid on
% box on
% hold on
% scatter(allK,allT2Profile,30,'Filled')
% scatter(allK,allT2ML,30,'Filled')
% %scatter(allK, T2am,30,'Filled')
% %scatter(allK, filteredT2ml,30,'Filled')
% %%scatter(allK, T2lm_long,30,'Filled')
% 
% 
% set(gca,'XScale','log')
% %set(gca,'YScale','log')
% 
% xlabel('K (m/s)')
% ylabel('T2 Estimator (s)')

%% 
figure(6)
grid on
box on 
hold on
scatter(allK, smallModeSum,30,'filled')
set(gca,'XScale','log')

%% 
figure(7)
grid on
box on
hold on
scatter(allT2ML, smallModeSum,30,'filled')

%%
% Small mode sum cutoff 0.04, T2ML cutoff 0.15
smallModeCutoff = 0.06;
T2MLmodeCutoff = 0.15;

KCutoff = 6.0*10^-5;


bestKestimate = zeros(length(allT2ML),1);
bestT2estimate = zeros(length(allT2ML),1);

for jj = 1:length(allT2ML)
    
    % Try removing long time data to see what happens
    longTimeCutoff = 0.15; %s
    [minVal, cutoffBin] = min(abs(T2linbins - longTimeCutoff));
    
    filteredT2Dist = allT2dist(:,1:cutoffBin);
    filteredSlice = filteredT2Dist(jj,:);
    filteredT2ml(jj,1) = 10.^(sum(T2logbins(1:cutoffBin).*filteredSlice)./sum(filteredSlice));
    
    
%     if smallModeSum(jj) > smallModeCutoff
%         bestT2estimate(jj,1) = filteredT2ml(jj,1);
%     else
%         bestT2estimate(jj,1) = allT2ML(jj,1);
%     end

    if allK(jj) < KCutoff
        bestT2estimate(jj,1) = filteredT2ml(jj,1);
    else
        bestT2estimate(jj,1) = allT2ML(jj,1);
    end
     
    
    
    
end

%%
figure(2)
grid on
box on
hold on
scatter(allK,allT2Profile,30,'Filled')
scatter(allK,allT2ML,30,'Filled')
%scatter(allK,bestT2estimate,30,'Filled')
%scatter(allK, T2am,30,'Filled')
%scatter(allK, filteredT2ml,30,'Filled')
%%scatter(allK, T2lm_long,30,'Filled')


set(gca,'XScale','log')
%set(gca,'YScale','log')

xlabel('K (m/s)')
ylabel('T2 Estimator (s)')

%%
KdiffFactor = estimateKdiffFactor(allK,allSDRKProfile,0);

figure()
grid on
box on
hold on
scatter(log10(allK),KdiffFactor)
xlim([min(log10(allK)) max(log10(allK))])
ylim([-1 30])
%set(gca,'XScale','log')
