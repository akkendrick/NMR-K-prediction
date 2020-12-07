% Try random forest classifer for data
clear


% Kansas
%sites = {'dpnmr_larned_east','dpnmr_larned_lwph','dpnmr_larned_west',...
%    'dpnmrA11','dpnmrA12','dpnmrC1S','dpnmrC1SE','dpnmrC1SW'};



% Washington
%sites = {'dpnmr_leque_east','dpnmr_leque_west'}
% 

% Wisconsin
% sites = {'Site1-WellG5','Site1-WellG6','Site2-WellPN1','Site2-WellPN2'};

 
 % all sites
 sites = {'Site1-WellG5','Site1-WellG6','Site2-WellPN1','Site2-WellPN2',...
   'dpnmr_larned_east','dpnmr_larned_lwph','dpnmr_larned_west',...
   'dpnmr_leque_east','dpnmr_leque_west','dpnmrA11','dpnmrA12',...
   'dpnmrC1S','dpnmrC1SE','dpnmrC1SW'};

 
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



%%
randOnes = ones(length(allT2ML),1);


%inputVars = [allT2ML T2hm T2am allPhi allGamma allEMI promRatio widthRatio allMinLoc smallModeSum largeModeSum];   
%inputNames = {'T2ML','T2hm','T2am','phi','gamma','EMI','prom','width','min','smallMode','largeMode'};

 %inputVars = [allT2ML T2hm T2am T2Seevers allPhi promRatio widthRatio allMinLoc smallModeSum largeModeSum allZ T2lm_short];
 
 inputVars = [T2am allPhi allMinLoc smallModeSum largeModeSum T2lm_short];
 %inputNames = {'T2ML','T2hm','T2am','T2Seevers','phi','prom','width','min','smallMode','largeMode' 'z'};
 
%inputVars = [ allT2ML T2hm T2am T2Seevers allPhi allMinLoc smallModeSum largeModeSum allZ];
%inputNames = {'T2ML','T2hm','T2am','T2Seevers','phi','min','smallMode','largeMode' 'z'};

inputNames = {'1', '2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'};



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

randInd = randi([1 225],75,1);

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
randInputVars = [randT2ML randT2hm randT2am randT2Seevers randPhi randProm randWidth randMinLoc randSmallModeSum randLargeModeSum randZ];
%randInputVars = [randT2ML randT2hm randT2am randT2Seevers randPhi randMinLoc randSmallModeSum randLargeModeSum randZ];


rng default

%tree = fitrtree(inputVars, allbProfile, 'CategoricalPredictors',10,'MinParentSize',30,'PredictorNames',inputNames);
tree = fitrtree(inputVars,allbProfile,'PredictorSelection','curvature','Surrogate','on');
%tree = fitrtree(inputVars,allbProfile,'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
%    struct('AcquisitionFunctionName','expected-improvement-plus'),'PredictorNames',inputNames)

%tree = fitrtree(inputVars,allK,'PredictorSelection','curvature','Surrogate','on','PredictorNames',inputNames);

%allSliceT2disc = [allSliceT2disc allPhi allZ];

%tree = fitrtree(allSliceT2disc,allbProfile,'PredictorSelection','curvature','Surrogate','on');
%tree = fitrtree(allT2dist,allbProfile,'PredictorSelection','curvature','Surrogate','on');


treeEnsemble = fitrensemble(inputVars, allbProfile, 'NumLearningCycles',300,'Method','LSBoost');

% treeEnsemble = regularize(treeEnsemble);

%%
% figure;
% semilogx(treeEnsemble.Regularization.Lambda,treeEnsemble.Regularization.ResubstitutionMSE, ...
%     'bx-','Markersize',10);
% line([1e-3 1e-3],[treeEnsemble.Regularization.ResubstitutionMSE(1) ...
%      treeEnsemble.Regularization.ResubstitutionMSE(1)],...
%     'Marker','x','Markersize',10,'Color','b');
% r0 = resubLoss(treeEnsemble);
% line([treeEnsemble.Regularization.Lambda(2) treeEnsemble.Regularization.Lambda(end)],...
%      [r0 r0],'Color','r','LineStyle','--');
% xlabel('Lambda');
% ylabel('Resubstitution MSE');
% annotation('textbox',[0.5 0.22 0.5 0.05],'String','unregularized ensemble', ...
%     'Color','r','FontSize',14,'LineStyle','none');
% %%
% ls = treeEnsemble;
% 
% rng(0,'Twister') % for reproducibility
% [mse,nlearn] = cvshrink(ls,'Lambda',ls.Regularization.Lambda,'KFold',5);
% 
% figure;
% semilogx(ls.Regularization.Lambda,ls.Regularization.ResubstitutionMSE, ...
%     'bx-','Markersize',10);
% hold on;
% semilogx(ls.Regularization.Lambda,mse,'ro-','Markersize',10);
% hold off;
% xlabel('Lambda');
% ylabel('Mean squared error');
% legend('resubstitution','cross-validation','Location','NW');
% line([1e-3 1e-3],[ls.Regularization.ResubstitutionMSE(1) ...
%      ls.Regularization.ResubstitutionMSE(1)],...
%     'Marker','x','Markersize',10,'Color','b','HandleVisibility','off');
% line([1e-3 1e-3],[mse(1) mse(1)],'Marker','o',...
%     'Markersize',10,'Color','r','LineStyle','--','HandleVisibility','off');
% %%
% figure;
% loglog(ls.Regularization.Lambda,sum(ls.Regularization.TrainedWeights>0,1));
% hold;
% 
% loglog(ls.Regularization.Lambda,nlearn,'r--');
% hold off;
% xlabel('Lambda');
% ylabel('Number of learners');
% legend('resubstitution','cross-validation','Location','NE');
% line([1e-3 1e-3],...
%     [sum(ls.Regularization.TrainedWeights(:,1)>0) ...
%     sum(ls.Regularization.TrainedWeights(:,1)>0)],...
%     'Marker','x','Markersize',10,'Color','b','HandleVisibility','off');
% line([1e-3 1e-3],[nlearn(1) nlearn(1)],'marker','o',...
%     'Markersize',10,'Color','r','LineStyle','--','HandleVisibility','off');
% %%
% jj = 1:length(ls.Regularization.Lambda);
% [jj;ls.Regularization.Lambda]

%%

%bPred = predict(tree,randInputVars);
%bPred = predict(tree, allSliceT2disc);
bPred = predict(tree, inputVars);
ensemblePred = predict(treeEnsemble,inputVars);
%ensemblePred = predict(treeEnsemble,allT2dist);

figure(1)
% scatter(randK, bPred)
% scatter(randK,ensemblePred)
scatter(allK, bPred)
scatter(allK,ensemblePred)
%
imp = predictorImportance(tree);

figure(3)
bar(imp);
title('Predictor Importance Estimates');
ylabel('Estimates');
xlabel('Predictors');
h = gca;
h.XTickLabel = tree.PredictorNames;
h.XTickLabelRotation = 45;
h.TickLabelInterpreter = 'none';

% Compute error in log space 
bPredLog = log10(bPred);
allbProfileLog = log10(allbProfile);


MAE = sum(abs(bPredLog - allbProfileLog))/length(bPredLog)

%%
cv = crossval(treeEnsemble,'KFold',5);
figure;
plot(kfoldLoss(cv,'Mode','Cumulative'));
xlabel('Number of trees');
ylabel('Cross-validated MSE');
%ylim([0.2,2])


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

