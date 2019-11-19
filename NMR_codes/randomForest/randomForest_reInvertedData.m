% Try random forest classifer for data
clear

% sites = {'dpnmr_larned_east','dpnmr_larned_lwph','dpnmr_larned_west',...
%    'dpnmr_leque_east','dpnmr_leque_west','dpnmrA11','dpnmrA12',...
%    'dpnmrC1S','dpnmrC1SE','dpnmrC1SW'};
% 
 sites = {'Site1-WellG5','Site1-WellG6','Site2-WellPN1','Site2-WellPN2'};

 
 % all sites
%  sites = {'Site1-WellG5','Site1-WellG6','Site2-WellPN1','Site2-WellPN2',...
%    'dpnmr_larned_east','dpnmr_larned_lwph','dpnmr_larned_west',...
%    'dpnmr_leque_east','dpnmr_leque_west','dpnmrA11','dpnmrA12',...
%    'dpnmrC1S','dpnmrC1SE','dpnmrC1SW'};

 
for kk = 1:length(sites)   
    
    smoothing = 1;
    [T2dist, T2linbins, nmrName] = loadInvertedNMRdata(sites{kk},smoothing);
    T2logbins = log10(T2linbins);
    %T2dist = T2dist(:,1:end-1);
    
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
%        
%        if ~isempty(gammaEMIdepth)
%            [minVal, index] = min(abs(gammaEMIdepth - z(jj)));
%           filtGammaIndices(jj,1) = index; 
%           filtGamma(jj,1) = gamma(index);        
%           filtGammaDepths(jj,1) = gammaEMIdepth(index);
%           filtEMI(jj,1) = EMI(index);
%        end
       
       %filtPhi(jj) = phi(index)
    end
    
    [sortedK{kk}, sortInd] = sort(K);
    sortedT2dist{kk} = filtT2dist(sortInd,:);
    sortedPhi{kk} = phi(sortInd);
    sortedT2ML{kk} = T2ML(sortInd);
    sortedZ{kk} = z(sortInd);
    sortedLogT2ML{kk} = logT2ML(sortInd);
    sortedLogK{kk} = logK(sortInd);

   % sortedGamma{kk} = filtGamma(sortInd);
    %sortedGammaDepths{kk} = filtGammaDepths(sortInd);
    %sortedEMI{kk} = filtEMI(sortInd);

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
    
    figure(20+kk)
    ribbon(sortedT2dist{kk}')
    
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
% Try plotting ribbon of all data (sorted)
[allSortedK, allSortInd] = sort(allK);
allSortedT2dist = allT2dist(allSortInd,:);

figure(25)
ribbon(allSortedT2dist')

disp('Nope')   
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

%     
%     T2am(jj,1) = sum(T2linbins.*sliceT2)./allPhi(jj); %Arithmetic Mean
%     T2lm(jj,1) = 10.^(sum(T2logbins.*sliceT2)./allPhi(jj));
%     %T2lm_filt(jj,1) = 10.^(sum(T2logbinsFilt.*sliceT2Filt)./allPhi(jj));
%     T2hm(jj,1) = allPhi(jj)./(sum(sliceT2./T2linbins)); 

    totSignal(jj,1) = sum(sliceT2);
    T2am(jj,1) = sum(T2linbins.*sliceT2)./totSignal(jj); %Arithmetic Mean
    T2lm(jj,1) = 10.^(sum(T2logbins.*sliceT2)./totSignal(jj));
    %T2lm_filt(jj,1) = 10.^(sum(T2logbinsFilt.*sliceT2Filt)./allPhi(jj));
    T2hm(jj,1) = totSignal(jj)./(sum(sliceT2./T2linbins)); 
    
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

    minInd = find(minLoc == 1);

    % Take minVals less than max
%    goodMinLoc = minInd(minInd < maxInd);
    goodMinLoc = minInd(minInd < 120);

    % Take the local min with the largest relax time
    goodMinLoc = max(goodMinLoc);
    

    if ~isempty(goodMinLoc)
        shortSummedSignal(jj) = sum(sliceT2(1:goodMinLoc));
        longSummedSignal(jj) = sum(sliceT2(goodMinLoc:end));
        allMinLoc(jj,1) = goodMinLoc;
        allMinT2(jj,1) = T2linbins(allMinLoc(jj,1));

        
        T2lm_short(jj,1) = 10.^(sum(T2logbins(1:goodMinLoc).*sliceT2(1:goodMinLoc))./shortSummedSignal(jj));

    else
        shortSummedSignal(jj) = 0;
        longSummedSignal(jj) = sum(sliceT2);
        allMinLoc(jj,1) = NaN;
        allMinT2(jj,1) = NaN;
        T2lm_short(jj,1) = NaN;
    end

    largeModeSum(jj,1) = longSummedSignal(jj);
    smallModeSum(jj,1) = shortSummedSignal(jj);
    
    
    promThresh = (max(sliceT2)-min(sliceT2))/10;
    
    % Compute max positions??
    maxLoc = islocalmax(sliceT2,'MinProminence',promThresh);
    maxInd = find(maxLoc == 1);
    allMaxInd(jj,1) = min(maxInd);
    allsubMaxT2(jj,1) = T2linbins(allMaxInd(jj,:));
    
    nBins = 200;
    binSize = length(sliceT2)/nBins;
    for kk = 1:nBins
        startInd = (kk-1)*binSize+1;
        endInd = kk*binSize;
        discretizedT2Amp(kk) = sum(sliceT2(startInd:endInd));
    end

    allDiscT2Amp(jj,:) = discretizedT2Amp;

    allT2cumsum(jj,:) = cumsum(sliceT2,'reverse');
    
    
end

SDR_b = @(K,m,n,phi,T2ML) K./((phi.^m).*(T2ML).^n);

m = 0;
n = 2;

%allbProfile = SDR_b(allK,m,n,allPhi,allT2ML);
allbProfile = SDR_b(allK,m,n,allPhi,allT2ML);
%%
sortedAllT2cumsum = allT2cumsum(allSortInd,:);

figure(26)
ribbon(sortedAllT2cumsum')

disp('Nope')   

%%
randOnes = ones(length(allT2ML),1);


%inputVars = [allT2ML T2hm T2am allPhi allGamma allEMI promRatio widthRatio allMinLoc smallModeSum largeModeSum];   
%inputNames = {'T2ML','T2hm','T2am','phi','gamma','EMI','prom','width','min','smallMode','largeMode'};

 inputVars = [T2lm T2hm T2am T2Seevers allPhi promRatio widthRatio allMinT2 allsubMaxT2 smallModeSum largeModeSum T2lm_short];
% inputNames = {'T2lm','T2hm','T2am','T2Seevers','phi','prom','width','min','subMaxT2','smallMode','largeMode', 'z', 'K'};
 
%inputVars = [ allT2ML T2hm T2am T2Seevers allPhi allMinLoc smallModeSum largeModeSum allZ];
%inputNames = {'T2ML','T2hm','T2am','T2Seevers','phi','min','smallMode','largeMode' 'z'};

%inputNames = {'1', '2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'};

% Try discretizing data
%inputVars = allT2cumsum;
%inputVars = allT2dist;
%inputNames = {'1', '2','3','4','5','6','7','8','9','10'};%,'11','12','13','14','15','16','17','18','19','20'};



figure(1)
hold on
grid on
box on

scatter(allK,allbProfile,30,'Filled')
%scatter(allK,allK,30,'Filled')
%scatter(allCorrK, allbProfile, 30,'filled')

set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'FontSize',14)
xlabel('K (m/s)')
ylabel('SDR b')
%ylabel('Random Forest K (m/s)')


% Take a subset of data for prediction
%randInd = randi([1 83],20,1);
%randInd = randi([1 225],225,1);

randInd = 1:1:225;
%randInd = 1:1:83;

randK = allK(randInd);
randT2ML = allT2ML(randInd);
randT2hm = T2hm(randInd);
randT2am = T2am(randInd);
randCorrT2ML = allCorrT2ML(randInd);


randPhi = allPhi(randInd);
% randGamma = allGamma(randInd);
% randEMI = allEMI(randInd);
randProm = promRatio(randInd);
randWidth = widthRatio(randInd);
randMinLoc = allMinLoc(randInd);
randSmallModeSum = smallModeSum(randInd);
randLargeModeSum = largeModeSum(randInd);
randT2Seevers = T2Seevers(randInd);
randZ = allZ(randInd);

randT2lm = T2lm(randInd);

%randInputVars = [randT2ML randT2hm randT2am randPhi randGamma randEMI randProm randWidth randMinLoc randSmallModeSum randLargeModeSum];
%randInputVars = [randT2ML randT2hm randT2am randT2Seevers randPhi randProm randWidth randMinLoc randSmallModeSum randLargeModeSum randZ];
%randInputVars = [randT2ML randT2hm randT2am randT2Seevers randPhi randMinLoc randSmallModeSum randLargeModeSum randZ];


rng default

%tree = fitrtree(inputVars, allbProfile, 'CategoricalPredictors',10,'MinParentSize',30,'PredictorNames',inputNames);

%tree = fitrtree(inputVars,allbProfile,'PredictorSelection','curvature','Surrogate','on','PredictorNames',inputNames);

%tree = fitrtree(inputVars,allbProfile,'PredictorSelection','curvature','Surrogate','on');
treeBag = TreeBagger(200,inputVars,allbProfile,'Method','regression','Surrogate','on',...
    'PredictorSelection','curvature','OOBPredictorImportance','on');
treeEnsemble = fitrensemble(inputVars, allbProfile);

%tree = fitrtree(inputVars,allbProfile,'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
%    struct('AcquisitionFunctionName','expected-improvement-plus'),'PredictorNames',inputNames)

%tree = fitrtree(inputVars,allK,'PredictorSelection','curvature','Surrogate','on'); %'PredictorNames',inputNames);

%allSliceT2disc = [allSliceT2disc allPhi allZ];

%tree = fitrtree(allSliceT2disc,allbProfile,'PredictorSelection','curvature','Surrogate','on');
%tree = fitrtree(allT2dist,allbProfile,'PredictorSelection','curvature','Surrogate','on');
%%

meanPred = predict(treeBag,inputVars);
quartilePred = quantilePredict(treeBag,inputVars,'Quantile',[0.25,0.5,0.75]);
ensemblePred = predict(treeEnsemble,inputVars);
%bPred = predict(tree, allSliceT2disc);
%bPred = predict(tree, allT2dist);

figure(1)
%scatter(randK, meanPred)
scatter(randK,ensemblePred)
%plot(randK,quartilePred,'*')
legend('Data','MeanPred','First Quartile','Median','Third Quartile')

%
%imp = predictorImportance(tree);
imp = treeBag.OOBPermutedPredictorDeltaError;

figure(3)
bar(imp);
title('Predictor Importance Estimates');
ylabel('Estimates');
xlabel('Predictors');
h = gca;
%h.XTickLabel = tree.PredictorNames;
h.XTickLabelRotation = 45;
h.TickLabelInterpreter = 'none';

% Compute error in log space 
bPredLog = log10(bPred);
allbProfileLog = log10(allbProfile);


MAE = sum(abs(bPredLog - allbProfileLog))/length(bPredLog)

%%
%optimalTree = fitrtree(inputVars,allbProfile,'MinLeafSize',27,'PredictorSelection','curvature','Surrogate','on','PredictorNames',inputNames);
% tree = fitrtree(allT2dist,allbProfile,'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
%    struct('AcquisitionFunctionName','expected-improvement-plus'))
% 
% optimalTree = fitrtree(allT2dist,allbProfile,'MinLeafSize',71,'PredictorSelection','curvature','Surrogate','on');
% 
% view(optimalTree,'Mode','graph')
% 
% 
% % New prediction
% bOptPred = predict(optimalTree,allT2dist);
% 
% figure(1)
% scatter(randK, bOptPred)
% 
% % Compute error in log space 
% bOptPredLog = log10(bOptPred);
% bPredLog = log10(bPred);
% allbProfileLog = log10(allbProfile);
% 
% 
% optMAE = sum(abs(bOptPredLog - allbProfileLog))/length(bOptPred)
% rmsError = sum(abs(bPredLog - allbProfileLog))/length(bOptPred)

%%
%view(tree,'Mode','graph')

%% Try to evaluate tree quality
% 
% resuberror = resubLoss(tree)
% 
% rng 'default'
% cvrtree = crossval(tree);
% cvloss = kfoldLoss(cvrtree)
% 
% % Don't seem to be super useful metrics
% 
% % Try trimming tree
% leafs = logspace(1,2,100);
% N = numel(leafs);
% err = zeros(N,1);
% for n=1:N
%     t = fitrtree(inputVars,allbProfile,'CrossVal','On',...
%         'MinLeafSize',leafs(n));
%     err(n) = kfoldLoss(t);
% end
% figure()
% plot(leafs,err);
% xlabel('Min Leaf Size');
% ylabel('cross-validated error');

