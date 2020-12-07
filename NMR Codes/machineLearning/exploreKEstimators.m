% Try random forest classifer for data
clear


% Kansas
%sites = {'dpnmr_larned_east','dpnmr_larned_lwph','dpnmr_larned_west',...
%    'dpnmrA11','dpnmrA12','dpnmrC1S','dpnmrC1SE','dpnmrC1SW'};



% Washington
%sites = {'dpnmr_leque_east','dpnmr_leque_west'}
% 

% Wisconsin
sites = {'Site1-WellG5','Site1-WellG6','Site2-WellPN1','Site2-WellPN2'};

 
 % all sites
% sites = {'Site1-WellG5','Site1-WellG6','Site2-WellPN1','Site2-WellPN2',...
 %  'dpnmr_larned_east','dpnmr_larned_lwph','dpnmr_larned_west',...
 %  'dpnmr_leque_east','dpnmr_leque_west','dpnmrA11','dpnmrA12',...
 %  'dpnmrC1S','dpnmrC1SE','dpnmrC1SW'};

 
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
    sortedSOE{kk} = SumEch(sortInd);
    
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
allSOE = vertcat(sortedSOE{:});

allT2linbins = vertcat(siteT2linbins{:});
allT2logbins = vertcat(siteT2logbins{:});
allCorrT2ML = 10.^(vertcat(corrlogT2ML{:}));
allCorrK = 10.^(vertcat(corrlogK{:}));
T2Seevers = (allT2ML.^(-1) - T2B^(-1)).^(-1);

   
%%
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

    
    T2am(jj,1) = sum(T2linbins.*sliceT2)./allPhi(jj); %Arithmetic Mean
    T2lm(jj,1) = 10.^(sum(T2logbins.*sliceT2)./allPhi(jj));
    T2hm(jj,1) = allPhi(jj)./(sum(sliceT2./T2linbins)); 
    
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

    minInd = find(minLoc == 1);

    % Take minVals less than max
    goodMinLoc = minInd(minInd < maxInd);
    
    % Take the local min with the largest relax time
    goodMinLoc = max(goodMinLoc);
    

    if ~isempty(goodMinLoc)
        shortSummedSignal(jj) = sum(sliceT2(1:goodMinLoc));
        longSummedSignal(jj) = sum(sliceT2(goodMinLoc:end));
        allMinLoc(jj,1) = goodMinLoc;
        
        T2lm_short(jj,1) = 10.^(sum(T2logbins(1:goodMinLoc).*sliceT2(1:goodMinLoc))./shortSummedSignal(jj));


    else
        shortSummedSignal(jj) = 0;
        longSummedSignal(jj) = sum(sliceT2);
        allMinLoc(jj,1) = NaN;
        T2lm_short(jj,1) = NaN;

    end

    largeModeSum(jj,1) = longSummedSignal(jj);
    smallModeSum(jj,1) = shortSummedSignal(jj);
    allT2cumsum(jj,:) = cumsum(sliceT2,'reverse');


end

SDR_b = @(K,m,n,phi,T2ML) K./((phi.^m).*(T2ML).^n);

m = 0;
n = 2;

%allbProfile = SDR_b(allK,m,n,allPhi,allT2ML);
allbProfile = SDR_b(allK,m,n,allPhi,allT2ML);
%%
% Plot T2 estimators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
box on
grid on
hold on

%scatter(allT2ML.^2, allK, 40,'Filled')
scatter(T2lm.^2, allK, 40, 'Filled')
scatter(T2am.^2, allK, 40, 'Filled')
scatter(T2hm.^2, allK, 40, 'Filled')

set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'FontSize',16)

ylim([10^-7,10^-2])
xlim([10^-4, 10^0])

ylabel('Hydraulic Conductivity (m/s)')
xlabel('T2.^2 (s)')

legend({'Log Mean','Arith Mean','Harmonic Mean'},'Location','northwest')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
box on
grid on
hold on

scatter(allSOE, allK, 40, 'Filled')

set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'FontSize',16)

xlabel('SOE')
ylabel('Hydraulic Conductivity (m/s)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)

subplot(1,2,1)
scatter(largeModeSum, allK, 40, 'Filled')
box on
grid on
set(gca,'YScale','log')
title('Large T2 Mode Sum')

subplot(1,2,2)
scatter(smallModeSum, allK, 40, 'Filled')
box on
grid on
set(gca,'YScale','log')
title('Small T2 Mode Sum')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)

testMinLoc = ~isnan(allMinLoc);

scatter(testMinLoc, allK, 40,'Filled')
box on
grid on
set(gca,'YScale','log')
title('Is T_2 dist bimodal?')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
T2B = 2.2293; % Avg pore water from Wisconsin
T2Seevers_vec = ((T2lm.^(-1) - T2B.^(-1)).^(-1));


% T2MLcorr = corrcoef(T2lm, allK)
% T2ML2corr = corrcoef(T2lm.^2, allK)
% T2amcorr = corrcoef(T2am.^2, allK)
% T2hmcorr = corrcoef(T2hm.^2, allK)
% SOEcorr = corrcoef(allSOE, allK)
% largeModeSumCorr = corrcoef(largeModeSum, allK)
% smallModeSumCorr = corrcoef(smallModeSum, allK)

corrMatrixGood = [T2lm.^2 T2am.^2 T2hm.^2 allSOE T2Seevers_vec.^2 allK];
corrMatrixBad = [largeModeSum smallModeSum promRatio widthRatio T2peakMax' T2peakMin' allK];
corrcoef(corrMatrixGood)
corrcoef(corrMatrixBad)


