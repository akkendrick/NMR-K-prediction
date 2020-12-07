% Estimate T2 parameters for all data 
clear
%close all


% sites = {'dpnmr_larned_east','dpnmr_larned_lwph','dpnmr_larned_west',...
%     'dpnmr_leque_east','dpnmr_leque_west','dpnmrA11','dpnmrA12',...
%     'dpnmrC1S','dpnmrC1SE','dpnmrC1SW','Site1-WellG5','Site1-WellG6',...
%     'Site2-WellPN1','Site2-WellPN2'};

sites = {'Site1-WellG5','Site1-WellG6','Site2-WellPN1','Site2-WellPN2'};

for kk = 1:length(sites)   
    
    [T2dist, T2logbins, nmrName] = loadRawNMRdata(sites{kk});
    
    gammaEMIdepth = [];
    [gammaEMIdepth, gamma, EMI] = loadGammaEMIData(sites{kk});
    
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
       
       if ~isempty(gammaEMIdepth)
           [minVal, index] = min(abs(gammaEMIdepth - z(jj)));
           filtGammaIndices(jj,1) = index; 
           filtGamma(jj,1) = gamma(index);        
           filtGammaDepths(jj,1) = gammaEMIdepth(index);
           filtEMI(jj,1) = EMI(index);
       end
       
       %filtPhi(jj) = phi(index)
    end

    [sortedK, sortInd] = sort(K);
    sortedT2dist = filtT2dist(sortInd,:);
    sortedPhi = phi(sortInd);
    sortedT2ML = T2ML(sortInd);
    z = z(sortInd)';
    [sortedLogT2ML, sortInd] = sort(logT2ML);
    sortedLogK = logK(sortInd,:);
    
    sortedGamma = filtGamma(sortInd);
    sortedGammaDepths = filtGammaDepths(sortInd);
    sortedEMI = filtEMI(sortInd);

    T2linbins = 10.^(T2logbins);
    
    
    % Compute various T2 parameters
    
    T2Seevers = (sortedT2ML.^(-1) - T2B^(-1)).^(-1);
    
    % Try applying long T2 correction
    T2MLlogCutoff = -0.9;
    
    smallLogK = sortedLogK(sortedLogT2ML < T2MLlogCutoff);
    largeLogK = sortedLogK(sortedLogT2ML > T2MLlogCutoff);
    smallLogT2ML = sortedLogT2ML(sortedLogT2ML < T2MLlogCutoff);
    largeLogT2ML = sortedLogT2ML(sortedLogT2ML > T2MLlogCutoff);
    
    corrFactor = 0.5;
    
    LogT2ML_corr = [smallLogT2ML; largeLogT2ML - corrFactor];
    T2ML_corr = 10.^LogT2ML_corr;
    
    T2am = 0;
    T2lm = 0;
    T2hm = 0;
    for jj = 1:length(z)
        sliceT2 = sortedT2dist(jj,:);
        T2am(jj,1) = sum(T2linbins.*sliceT2)./sortedPhi(jj); %Arithmetic Mean
        T2lm(jj,1) = 10.^(sum(T2logbins.*sliceT2)./sortedPhi(jj));
        T2hm(jj,1) = sortedPhi(jj)./(sum(sliceT2./T2linbins)); 
    end    

    promRatio = 0;
    widthRatio = 0;
    largeModeSum = 0;
    smallModeSum = 0;
    % Compute T2 distribution statistics
    for jj = 1:length(sortedT2ML)
        sliceT2 = sortedT2dist(jj,:);
       
        [pks{jj}, locs{jj},widths{jj},proms{jj}] = findpeaks(sliceT2);

        if length(pks{jj}) > 1
            T2peakAvg(jj) =  mean(T2linbins(locs{jj}));
            currentProm = proms{jj};
            currentWidth = widths{jj};
            currentLocs = locs{jj};

            promRatio(jj) = currentProm(1)/currentProm(2);
            widthRatio(jj) = currentWidth(1)/currentWidth(2);

            T2peakMax(jj) = T2linbins(currentLocs(2));
            T2peakMin(jj) = T2linbins(currentLocs(1));
        else
            T2peakAvg(jj) = T2linbins(locs{jj});
            promRatio(jj) = NaN;
            widthRatio(jj) = NaN;
            T2peakMax(jj) = T2linbins(locs{jj});
            T2peakMin(jj) = NaN;
        end
        
        %largeModeSum(jj,1) = sum(sliceT2(T2linbins >= 0.054));
        %smallModeSum(jj,1) = sum(sliceT2(T2linbins < 0.054));
        
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
        else
            shortSummedSignal(jj) = 0;
            longSummedSignal(jj) = sum(sliceT2);
        end
        
        largeModeSum(jj,1) = longSummedSignal(jj);
        smallModeSum(jj,1) = shortSummedSignal(jj);
        
    end
    
    T2_lowerLim = 0;
    T2_upperLim = 0;
    for jj = 1:length(sortedT2ML)
        if (widthRatio(jj) > mean(widthRatio,'omitNan'))...
                && (promRatio(jj) > mean(promRatio,'omitNaN'))...
                && (smallModeSum(jj) > mean(smallModeSum))...
                && (largeModeSum(jj) < mean(largeModeSum))            
            
            T2_lowerLim(jj,1) = T2peakMin(jj);
            T2_upperLim(jj,1) = T2am(jj);
        else
            T2_lowerLim(jj,1) = T2am(jj);
            T2_upperLim(jj,1) = T2am(jj);
        end
    end
 
    SDR_b = @(K,m,n,phi,T2ML) K./((phi.^m).*(T2ML).^n);

    m = 0;
    n = 1;
   
    bProfileT2ML{kk} = SDR_b(sortedK,m,n,sortedPhi,T2lm);
    bProfileT2lower{kk} = SDR_b(sortedK,m,n,sortedPhi,T2_lowerLim);
    bProfileT2MLCorr{kk} = SDR_b(sortedK,m,n,sortedPhi,T2Seevers);
    
    promRatio_cells{kk} = promRatio';
    widthRatio_cells{kk} = widthRatio';
    smallModeSum_cells{kk} = smallModeSum;
    largeModeSum_cells{kk} = largeModeSum;
    
    K_cells{kk} = sortedK;
    T2ML_cells{kk} = sortedT2ML;
    T2am_cells{kk} = T2am;
    T2lm_cells{kk} = T2lm;
    T2hm_cells{kk} = T2hm;
    phi_cells{kk} = sortedPhi;
    
    T2opt_cells{kk} = T2_lowerLim;
    T2ML_corr_cells{kk} = T2ML_corr;
    
    gamma_cells{kk} = sortedGamma;
    gammaDepths_cells{kk} = sortedGammaDepths;
    EMI_cells{kk} = sortedEMI;
    
      
end

%%

allbProfileT2ML = vertcat(bProfileT2ML{:});
allbProfileT2lower = vertcat(bProfileT2lower{:});
allbProfileT2ML_corr = vertcat(bProfileT2MLCorr{:});

allPromRatio = vertcat(promRatio_cells{:});
allWidthRatio = vertcat(widthRatio_cells{:});
allSmallModeSum = vertcat(smallModeSum_cells{:});
allLargeModeSum = vertcat(largeModeSum_cells{:});

allK = vertcat(K_cells{:});
allT2ML = vertcat(T2ML_cells{:});
allT2am = vertcat(T2am_cells{:});
allT2lm = vertcat(T2lm_cells{:});
allT2hm = vertcat(T2hm_cells{:});
allPhi = vertcat(phi_cells{:});
allT2opt = vertcat(T2opt_cells{:});
allT2ML_corr = vertcat(T2ML_corr_cells{:});

allGamma = vertcat(gamma_cells{:});
allGammaDepths = vertcat(gammaDepths_cells{:});
allEMI = vertcat(EMI_cells{:});

figure(2)
box on
grid on
hold on
Krange = [10^-6 10^-3];

scatter(allK, allbProfileT2ML,30,'Filled')
%scatter(allK, allbProfileT2lower,30,'Filled')
%scatter(allK, allbProfileT2ML_corr,30,'Filled')

xlabel('Hydraulic Conductivity (m/s)')
ylabel('Ideal SDR b parameter')

set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'FontSize',14)

xlim([10^-7 10^-2])
ylim([10^-5 10^1])

%%
% Try fitting hm and b data

allT2hm = allT2lm;

X = [ones(length(allT2hm),1) log10(allT2hm)];
fitParams = X\(log10(allbProfileT2ML));

T2hm_range = logspace(-2.4,-0.4,200);
T2hm_range = log10(T2hm_range);
bestFitb = fitParams(1) + (T2hm_range.*fitParams(2));

T2hm_linRange = 10.^T2hm_range;
bestFitb_lin = 10.^bestFitb;

otherBestFitb = 10^(-2.6985).*(T2hm_linRange.^(-0.6722));

figure(3)
box on
grid on
hold on

%scatter(allT2hm, allbProfileT2ML_corr,30,'Filled')
scatter(allT2hm, allbProfileT2ML,30,'Filled')
scatter(T2hm_linRange, bestFitb_lin,30,'Filled')
%scatter(T2hm_linRange, otherBestFitb, 30,'Filled')

set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'FontSize',14)

xlabel('T2hm (s)')
ylabel('Ideal SDR b parameter')




%%

% figure(4)
% box on
% grid on
% hold on
% 
% scatter(allPhi, allbProfileT2ML,30,'Filled')
% 
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% set(gca,'FontSize',14)
% 
% xlabel('Phi')
% ylabel('Ideal SDR b parameter')


% figure(5)
% box on
% grid on
% hold on
% 
% scatter(allWidthRatio, allbProfileT2ML, 30,'Filled')
% 
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% set(gca,'FontSize',14)
% 
% xlabel('Prom Ratio')
% ylabel('Ideal SDR b parameter')
%  

figure(6)
box on
grid on
hold on

scatter(allLargeModeSum,allbProfileT2ML,30,'Filled')

%set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'FontSize',14)

xlabel('Small Mode Sum')
ylabel('Ideal SDR b parameter')

figure(7)
box on
grid on
hold on

scatter(allGamma, allbProfileT2ML,30,'Filled')
set(gca,'yscale','log')
set(gca,'FontSize',14)

xlabel('Gamma')
ylabel('Ideal SDR b parameter')
%xlim([0,50])
