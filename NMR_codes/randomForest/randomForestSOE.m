% SOE Random Forest Classifier 

% Try random forest classifer for data
clear


%%Kansas
% sites = {'dpnmr_larned_east','dpnmr_larned_lwph','dpnmr_larned_west',...
%    'dpnmrA11','dpnmrA12','dpnmrC1S','dpnmrC1SE','dpnmrC1SW'};
% 
% 

% Washington
%sites = {'dpnmr_leque_east','dpnmr_leque_west'}
% 

% Wisconsin
% sites = {'Site1-WellG5','Site1-WellG6','Site2-WellPN1','Site2-WellPN2'};

 
%  % all sites
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
    sortedSOE{kk} = SumEch(sortInd);

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
allSOE = vertcat(sortedSOE{:});
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

m = 0;

SDR_b = @(K,m,n,phi,T2ML) K./((phi.^m).*(T2ML).^n);
SOE_b = @(K,n,SOE) K ./ (SOE.^n);

allSDRbProfile = SDR_b(allK,m,2,allPhi,allT2ML);
allSOEbProfile = SOE_b(allK,1,allSOE);

%%
inputVars = [T2am T2lm T2hm allPhi allSOE];
inputNames = {'1', '2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'};

figure(1)
hold on
grid on
box on

scatter(allK,allSOEbProfile,30,'Filled')
%scatter(allK,allK,30,'Filled')
%scatter(allCorrK, allbProfile, 30,'filled')

set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'FontSize',14)
xlabel('K (m/s)')
ylabel('SOE b')
%ylabel('Random Forest K (m/s)')

rng default

%tree = fitrtree(inputVars, allbProfile, 'CategoricalPredictors',10,'MinParentSize',30,'PredictorNames',inputNames);
tree = fitrtree(inputVars,allSOEbProfile,'PredictorSelection','curvature','Surrogate','on');
%tree = fitrtree(inputVars,allbProfile,'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
%    struct('AcquisitionFunctionName','expected-improvement-plus'),'PredictorNames',inputNames)

%tree = fitrtree(inputVars,allK,'PredictorSelection','curvature','Surrogate','on','PredictorNames',inputNames);

%allSliceT2disc = [allSliceT2disc allPhi allZ];

%tree = fitrtree(allSliceT2disc,allbProfile,'PredictorSelection','curvature','Surrogate','on');
%tree = fitrtree(allT2dist,allbProfile,'PredictorSelection','curvature','Surrogate','on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Pull random sample of data for ensemble fit
NsubSamp = 150;

randInd = randperm(225,NsubSamp);
inputVarsRand = inputVars(randInd,:);
allSOEbProfileRand = allSOEbProfile(randInd);

%treeEnsemble = fitrensemble(inputVars, allSOEbProfile, 'NumLearningCycles',300,'Method','LSBoost');
treeEnsemble = fitrensemble(inputVarsRand, allSOEbProfileRand, 'NumLearningCycles',300,'Method','LSBoost');


%%
bPred = predict(tree, inputVars);
ensemblePred = predict(treeEnsemble,inputVars);

figure(1)
%scatter(allK, bPred)
scatter(allK,ensemblePred)

xlim([10^-7, 10^-2])
ylim([10^-6, 10^1])

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
allbProfileLog = log10(allSOEbProfile);

ensemblePredLog = log10(ensemblePred);

MAE = sum(abs(bPredLog - allbProfileLog))/length(bPredLog)
MAE_ensemble = sum(abs(bPredLog - ensemblePredLog))/length(bPredLog)


%%
cv = crossval(treeEnsemble,'KFold',5);
figure;
plot(kfoldLoss(cv,'Mode','Cumulative'));
xlabel('Number of trees');
ylabel('Cross-validated MSE');
%ylim([0.2,2])


% %%
% 
% inputNames = {'T2am','T2lm','T2hm','allPhi','SOE'}
% 
% %optimalTree = fitrtree(inputVars,allSOEbProfile,'MinLeafSize',27,'PredictorSelection','curvature','Surrogate','on','PredictorNames',inputNames);
% 
% 
% %optimalTree = fitrtree(inputVars,allSOEbProfile,'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
% %   struct('AcquisitionFunctionName','expected-improvement-plus'))
% 
% optimalTree = fitrtree(allT2dist,allSOEbProfile,'MinLeafSize',71,'PredictorSelection','curvature','Surrogate','on');
% 
% view(optimalTree,'Mode','graph')
% 
% %%
% % New prediction
% bOptPred = predict(optimalTree,inputVars);
% 
% figure(1)
% scatter(allK, bOptPred)
% 
% % Compute error in log space 
% bOptPredLog = log10(bOptPred);
% bPredLog = log10(bPred);
% allbProfileLog = log10(allbProfile);
% 
% 
% optMAE = sum(abs(bOptPredLog - allbProfileLog))/length(bOptPred)
% rmsError = sum(abs(bPredLog - allbProfileLog))/length(bOptPred)

