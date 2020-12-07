% % Wisconsin
sites = {"Site1-WellG5","Site1-WellG6","Site2-WellPN1","Site2-WellPN2"};

% SDR_b = [0.0114 0.0065 0.0212 0.0132];
% SDR_m = [1 1 1 1];
% SDR_n = [2 2 2 2];

SDR_b = [0.0040 0.0025 0.0075 0.0047];
SDR_m = [0 0 0 0];
SDR_n = [2 2 2 2];


SDR_K = @(b,m,n,phi,T2ML) b.*(phi.^m).*(T2ML).^n;

% SDR_b_arith = [0.0072 0.0051 0.0190 0.0103];
% SDR_b_harm = [0.0471 0.0176 0.0319 0.0376];
% SDR_b_geom = [0.0112 0.0066 0.0213 0.0132];

SDR_b_arith = [0.0040 0.0025 0.0075 0.0047];
SDR_b_harm = [0.0040 0.0025 0.0075 0.0047];
SDR_b_geom = [0.0040 0.0025 0.0075 0.0047];

% SDR_b_arith = [5.8320 7.8231 23.6286 14.1389];
% SDR_b_harm = [0.0299 0.0104 0.0191 0.0193];
% SDR_b_geom = [0.4159 0.2053 0.03234 0.5051];


figureson = 0; 

for kk = 1:length(sites)   
    
    [T2dist, T2logbins, nmrName] = loadRawNMRdata(sites{kk});
    
   % gammaEMIdepth = [];
   % [gammaEMIdepth, gamma, EMI] = loadGammaEMIData(sites{kk});
    
    [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(nmrName);

    T2depths = T2dist(:,1);
    T2data = T2dist(:,2:end);
    T2linbins = 10.^T2logbins;
    
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
    currentSDRb_arith{kk} = ones(length(phi),1).*SDR_b_arith(kk);
    currentSDRb_harm{kk} = ones(length(phi),1).*SDR_b_harm(kk);
    currentSDRb_geom{kk} = ones(length(phi),1).*SDR_b_geom(kk);

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
allSDRb_arith = vertcat(currentSDRb_arith{:});
allSDRb_harm = vertcat(currentSDRb_harm{:});
allSDRb_geom = vertcat(currentSDRb_geom{:});

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
allSDRb_arith = allSDRb_arith(allSortInd);
allSDRb_harm = allSDRb_harm(allSortInd);
allSDRb_geom = allSDRb_geom(allSortInd);

%% %Compute different T2ML metric
for jj = 1:length(allT2ML)
    
    sliceT2 = allT2dist(jj,:);
    
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
    
%      T2am(jj,1) = sum(T2linbins.*sliceT2)./allPhi(jj); %Arithmetic Mean
%     T2hm(jj,1) = allPhi(jj)./(sum(sliceT2./T2linbins)); 
    T2hm_quant = 0;
    gmcount = 0;
    gmprod = 1;
    
    sliceT2_filt = ones(1,length(sliceT2));
    for kk = 1:1:length(sliceT2)
        %Assuming T2bins in linear space
         T2hm_quant = T2hm_quant + sliceT2(kk)/T2linbins(kk);
         
         T2_prod(1,kk) = sliceT2(1,kk)*T2linbins(kk); 
         T2_prod_log(1,kk) = sliceT2(1,kk)*T2linbins(kk); 
         %Is sliceT2 really the quantity we want here?
         if sliceT2(1,kk) ~= 0
             T2_filt(1,kk) = sliceT2(1,kk);
             gmcount = gmcount + 1;
         else
            T2_prod(1,kk) = 1;
            T2_prod_log(1,kk) = 1; 
         end
    end
    
    
    % Need to fix T2gm problems
    % Update 7/25
    % I thought I fixed it this time, but I haven't. What is still wrong
    % here or will it always not work?
    
    T2gm_quant = cumprod(T2_prod,2);
    
    %T2am(jj,1) = sum(sliceT2)/length(sliceT2);
    T2am_weight(jj,1) = sum(T2linbins.*sliceT2)./allPhi(jj);

%     T2gm(jj,1) = T2gm_quant(end)^(1/gmcount);
%     T2gmlog(jj,1) = exp(sum(log(sliceT2_filt))*(1/gmcount));
    
    T2gm_weight(jj,1) = T2gm_quant(end)^(1/gmcount);
    T2gm_weightlog(jj,1) = exp(sum(log(T2_prod_log))*(1/gmcount));
    
    T2hm_weight(jj,1) = sum(sliceT2)/T2hm_quant;
    T2lm(jj,1) = 10.^(sum(T2logbins.*sliceT2)./allPhi(jj));

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
            
%       This arithmetic averaging gives the same result  
%         arithT2ML(jj,1) = (longSummedSignal(jj)/(longSummedSignal(jj) + shortSummedSignal(jj)))...
%             *T2lm_long(jj,1) + (shortSummedSignal(jj)/(longSummedSignal(jj) + shortSummedSignal(jj)))*T2lm_short(jj,1);
   
%         arithT2ML(jj,1) = (T2lm_short(jj)*shortSummedSignal(jj))/(allPhi(jj))+...
%             (T2lm_long(jj,1)*longSummedSignal(jj))/allPhi(jj);
%         
%         harmT2ML(jj,1) = (((longSummedSignal(jj)/T2lm_long(jj,1)+...
%             shortSummedSignal(jj)/T2lm_short(jj,1)))/allPhi(jj))^(-1);
%         
%         geomT2ML(jj,1) = ((T2lm_long(jj,1)^longSummedSignal(jj))*...
%             T2lm_short(jj,1)^shortSummedSignal(jj))^(1/allPhi(jj));
      

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
        T2lm_long(jj,1) = allT2ML(jj);
        
        T2_bimodal(jj) = 0;
        
        [maxVal, indLongMax(jj)] = max(longSignal);
        longT2Max(jj,1) = T2linbins(indLongMax(jj));
        shortT2Max(jj,1) = 0;
        indShortMax(jj,1) = 1;
        
%         arithT2ML(jj,1) = (T2lm_long(jj,1)*longSummedSignal(jj))/allPhi(jj);
%         
%         harmT2ML(jj,1) = ((longSummedSignal(jj)/T2lm_long(jj,1))/...
%             allPhi(jj))^(-1);
%         
%         geomT2ML(jj,1) = ((T2lm_long(jj,1)^longSummedSignal(jj)))^(1/allPhi(jj));
        
        

        
    end 
        if figureson == 1
            allK(jj)

            %Checked on 6/6 T2ML and T2lm are the same
            [val, index] = min(abs(allT2ML(jj) - T2linbins));
            [val, index2] = min(abs(arithT2ML(jj) - T2linbins));
            [val, index3] = min(abs(harmT2ML(jj) - T2linbins));
            [val, index4] = min(abs(geomT2ML(jj) - T2linbins));

            plot(T2linbins, sliceT2, 'LineWidth',2)

            grid on
            box on
            hold on

            plot(T2linbins, longSignals(1:100,jj),'--', 'LineWidth',2)
            plot(T2linbins, shortSignals(1:100,jj),':', 'LineWidth',2)

            plot(T2linbins(index), sliceT2(index), 'k*','MarkerSize',10)
            plot(T2linbins(index2), sliceT2(index2), 'r*','MarkerSize',10)
            plot(T2linbins(index3), sliceT2(index3), 'g*','MarkerSize',10)
            plot(T2linbins(index4), sliceT2(index4), 'b*','MarkerSize',10)

            smallwtstr = num2str(shortSummedSignal(jj));
            largewtstr = num2str(longSummedSignal(jj));
            
            title(['Short weight=',smallwtstr, ' Long weight=',largewtstr, ' phi=',num2str(allPhi(jj))])
    %         plot(T2linbins(index), sliceT2(index), 'k*','MarkerSize',10)
%             plot(T2linbins(indLongMax(jj,1)), sliceT2(indLongMax(jj,1)), 'r*','MarkerSize',10)
%             plot(T2linbins(indShortMax(jj,1)), sliceT2(indShortMax(jj,1)), 'g*','MarkerSize',10)


            hold off
            set(gca,'XScale','log')
            
            legendString = ['Long Signal','Short Signal','T2ML','Arith T2ML'...
                ,'Harmonic T2ML','Geometric T2ML'];
            %legendString = strcat('K= ', sprintf('%d',allK(jj)),'m/s');
            %fileString = strcat('K= ', sprintf('%d',allK(jj)), allSites(jj),'.png');

            fileString = strcat('z = ',num2str(allZ(jj)),'K=',sprintf('%d',allK(jj)), allSites(jj),'.png');

            legend('T2dist','Long Signal','Short Signal','T2ML','Arith T2ML'...
                ,'Harmonic T2ML','Geometric T2ML','Location','northwest')
            print('-dpng','-r300',fileString)
        end
end

for jj = 1:length(allT2ML)
    optSDR_K(jj) = SDR_K(allSDRb(jj), allSDRm(jj), allSDRn(jj),...
    allPhi(jj), allT2ML(jj));

    arithSDR_K(jj) = SDR_K(allSDRb_arith(jj), allSDRm(jj), allSDRn(jj),...
    allPhi(jj), T2am_weight(jj));

    harmSDR_K(jj) = SDR_K(allSDRb_harm(jj), allSDRm(jj), allSDRn(jj),...
    allPhi(jj), T2hm_weight(jj));

    geomSDR_K(jj) = SDR_K(allSDRb_geom(jj), allSDRm(jj), allSDRn(jj),...
    allPhi(jj), T2gm_weight(jj));
end

optKdiff = median(estimateKdiffFactor(allK,optSDR_K',1))
arithKdiff = median(estimateKdiffFactor(allK,arithSDR_K',1))
harmKdiff = median(estimateKdiffFactor(allK,harmSDR_K',1))
geomKdiff= median(estimateKdiffFactor(allK,geomSDR_K',1))

% [optSign, optKdiff] = estimateKdiffFactor_withSign(allK,optSDR_K,1);
% [arithSign, arithKdiff] = estimateKdiffFactor_withSign(allK,arithSDR_K,1);
% [harmSign, harmKdiff] = estimateKdiffFactor_withSign(allK,harmSDR_K,1);
% [geomSign, geomKdiff] = estimateKdiffFactor_withSign(allK,geomSDR_K,1);


k_estimates = [optSDR_K; arithSDR_K; harmSDR_K; geomSDR_K];
%k_estimates = [optSDR_K; arithSDR_K; harmSDR_K];

k_sym = ['+','.','.','.'];
k_size = [22,20,20,20];
k_estimates = k_estimates';
k_names = {'Regular SDR T2ML','Arithmetic T2 Mean', 'Harmonic T2 Mean' 'Geometric T2 Mean'};
plotKestKdpp_sym(allK,k_estimates,k_estimates,k_names,k_sym,k_size)