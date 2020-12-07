% Estimate porosity weighted T2 value

%clear
%close all

clear 
% sites = {'dpnmr_larned_east','dpnmr_larned_lwph','dpnmr_larned_west',...
%     'dpnmr_leque_east','dpnmr_leque_west','dpnmrA11','dpnmrA12',...
%     'dpnmrC1S','dpnmrC1SE','dpnmrC1SW'};
% 

%sites = {'Site1-WellG5','Site1-WellG6','Site2-WellPN1','Site2-WellPN2'};

sites = {'Site1-WellG5'};
legendNames = [{'Site1-WellG6'}];

waterTable = [2.1248,2.0469,5.0285,4.7476]; % rel ground surface
depthOffsets = [0.95,0.75,0.75,0.75];

for kk = 1:length(sites)
    %close all
    baseName = sites{kk}
%     
    rawBaseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/USGS Data/';
    baseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/NMR-K-prediction/Data/Aggregated_Data/';
    
      
    %rawBaseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/Kansas_Wash_Data/';
    %baseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/NMR-K-prediction/Data/Aggregated_Data/';
    
    
    
%     rawBaseDir = 'I:\My Drive\USGS Project\USGS Data\';
%     baseDir = 'I:\My Drive\USGS Project\NMR-K-prediction/Data/Aggregated_Data/';
    
    if strcmp(baseName,'Site1-WellG5')
        site = 'Site1-WellG5';
        name = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1';
        nmrName = name;
    elseif  strcmp(baseName,'Site1-WellG5above')
        site = 'Site1-WellG5';
        name = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1';
        nmrName = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1_above';
    elseif  strcmp(baseName,'Site1-WellG5below')
        site = 'Site1-WellG5';
        name = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1';
        nmrName = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1_below';
    elseif strcmp(baseName,'Site1-WellG6')
        site = 'Site1-WellG6'
        name = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';
        nmrName = name;
    elseif strcmp(baseName,'Site1-WellG6above')
        site = 'Site1-WellG6';
        name = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';
        nmrName = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1_above'
    elseif strcmp(baseName,'Site1-WellG6below')
        site = 'Site1-WellG6';
        name = 'G6_W2_ts5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';
        nmrName = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1_below'
    elseif strcmp(baseName,'Site2-WellPN1')
        site = 'Site2-WellPN1';
        name = 'Pl_W1_Tr5_20x_MPp75aLS_F1n2_wRIN_wRFI_Reg50_Va1';
        nmrName = name;
    elseif strcmp(baseName,'Site2-WellPN2')
        site = 'Site2-WellPN2';
        name = 'W2_Tr5_20x_MPp75aLS_Reg50_wRIN_wRFI_Va1';
        nmrName = name;
    else
        site = baseName;
        nmrName = baseName;
        name = baseName;
    end

    [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(nmrName);

    in1 = [rawBaseDir site '/' name '/' name '_T2_dist' '.txt']; 
    in2 = [rawBaseDir site '/' name '/' name '_T2_bins_log10s' '.txt']; 
    in3 = [rawBaseDir site '/' name '/' name '_1Dvectors' '.txt'];

%     in1 = [rawBaseDir site '/' name '_T2_dist' '.txt']; 
%     in2 = [rawBaseDir site '/' name '_T2_bins_log10s' '.txt']; 
%     
    T2dist = load(in1);
    %T2logbins = load(in2);
    T2logbins = load(in2);
    
    %dataVectors = dlmread(in3,'',1,0);
    
    %phiVector = dataVectors(:,2);
    
    T2dist(:,1) = T2dist(:,1) - depthOffsets(2); 
       
    T2depths = T2dist(:,1);
    T2data = T2dist(:,2:end);
    
    T2axis = 10.^T2logbins;
    T2axis = T2axis * 10^3;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find T2 distributions where we have DPP data
    
    for jj = 1:length(z)
       [minVal, index] = min(abs(T2depths - z(jj)));
       filtIndices(jj) = index; 
       filtT2dist(jj,:) = T2data(index,:);        
       filtT2depths(jj) = T2depths(index);
       %filtPhi(jj) = phiVector(index);
    end
        
    %filtPhi = phi;
    [sortedK, sortInd] = sort(K);
    sortedT2dist = filtT2dist(sortInd,:);
    %sortedPhi = filtPhi(sortInd)';
    sortedPhi = phi(sortInd);
    T2ML = T2ML(sortInd)';
    z = z(sortInd)';
    
    figure(20)
    ribbon(filtT2dist')
    title(baseName)
    
    %%
%     % Quickly look at porosity estimator
%     sliceT2 = filtT2dist(5,:);
%     slicePhi = filtPhi(5);
%     
%     plot(T2logbins, sliceT2)
%     %set(gca,'XScale','log')
%     
    %%
    % Now estimate porosity weighted T2 
    T2linbins = 10.^(T2logbins);

    for jj = 1:length(z)
        sliceT2 = sortedT2dist(jj,:);
        T2W(jj,1) = sum(T2linbins.*sliceT2)./sortedPhi(jj);
    end    
  
    figure(3)
    hold on
    
    plot(T2W)
    plot(T2ML)
    %%
    
    T2estimator = T2W;
    
    Kcm = sortedK *10^2; % m/s -> cm/s
    
    voidRatio = sortedPhi ./ (1-sortedPhi);
    refVoidRatio = 1;
    
    beta = 4;
    refK = Kcm.*(voidRatio./refVoidRatio).^(-beta);
    
    logSpecSurf = (log10(refK)+5)./(-2);
    specSurf = 10.^(logSpecSurf); % in m^2/g
           
    rho_m = 2650; % kg/m^3 for quartz, seems representative
    
    meanPoreSize = 2*(voidRatio)./(specSurf.*rho_m); %in mm
%     meanPoreSizeMM = meanPoreSize.*10^3;
    
    figure(1) 
    grid on
    hold on
    
    scatter(meanPoreSize,Kcm,40,'Filled')
    
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    
    xlabel('Estimated mean pore size (mm)')
    ylabel('Hydraulic Conductivity (cm/s)')
    
    box on
    
    set(gca,'FontSize',12)
    
    figure(2)
    grid on
    hold on
    
    scatter(meanPoreSize, T2estimator.^2,40,'Filled')
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    
    xlabel('Estimated mean pore size (mm)')
    ylabel('T_{2ML}^2 (s)')
    
    box on
    
    set(gca,'FontSize',12)
    
      
    % Attempt to estimate pore radius from Godefroy
    D = @(temp) (1.0413+0.039828*temp+0.00040318*temp^2)*10^-9;
    T2B = 2.2;
    alpha = 3;
    %bParam = (T2ML.*T2B.*alpha)./(T2B - T2ML);
    
    SeeversT2 = (T2ML.^(-1) - T2B^(-1)).^(-1);
    SeeversT2_take2 = (T2ML.*T2B)./(T2B-T2ML);
    
    bParam = (T2estimator.*T2B.*alpha)./(T2B - T2estimator);
    %bParam = (SeeversT2.*T2B.*alpha)./(T2B-SeeversT2);
    
    
    
    rho = 44*10^-6; %m/s
    rhoMod = 10 ;
    rho = rho * rhoMod;
    
    Dtemp = D(15);
    
    poreSize = (sqrt(Dtemp.*(2.*bParam.*rho^2+Dtemp)) - Dtemp)./rho;
    
    figure(2)
    scatter(poreSize*10^3, T2estimator.^2,30,'Filled')

    figure(1) 
    scatter(poreSize*10^3, Kcm, 30, 'Filled') 
    
%     % Attempt to estimate an r from rho S/V
%     rEst = (3*T2estimator*rho);
%     dEst = 2*rEst;
%     
%     figure(2)
%     scatter(dEst*10^3, T2estimator.^2,30,'Filled')

%%
% Attempt calculating pore size for each T2 mode

 
for jj = 1:length(T2estimator)
   sliceT2 = sortedT2dist(jj,:);
   bParam = (T2linbins.*T2B.*alpha)./(T2B - T2linbins);
   poreSizeDist = ((sqrt(Dtemp.*(2.*bParam.*rho^2+Dtemp)) - Dtemp)./rho)*10^3;    

   T2MLbParam = (T2ML(jj).*T2B.*alpha)./(T2B - T2ML(jj));
   T2ML_poreSize(jj) = ((sqrt(Dtemp.*(2.*T2MLbParam.*rho^2+Dtemp)) - Dtemp)./rho)*10^3; 
      
   T2WbParam = (T2W(jj).*T2B.*alpha)./(T2B - T2W(jj));
   T2W_poreSize(jj) = ((sqrt(Dtemp.*(2.*T2WbParam.*rho^2+Dtemp)) - Dtemp)./rho)*10^3; 
   
   
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
   
   T2peakAvgbParam = (T2peakAvg(jj).*T2B.*alpha)./(T2B - T2peakAvg(jj));
   T2peakAvg_poreSize(jj) = ((sqrt(Dtemp.*(2.*T2peakAvgbParam.*rho^2+Dtemp)) - Dtemp)./rho)*10^3; 
   
   T2peakMaxbParam = (T2peakMax(jj).*T2B.*alpha)./(T2B - T2peakAvg(jj));
   T2peakMax_poreSize(jj) = ((sqrt(Dtemp.*(2.*T2peakMaxbParam.*rho^2+Dtemp)) - Dtemp)./rho)*10^3; 
   
   T2peakMinbParam = (T2peakMin(jj).*T2B.*alpha)./(T2B - T2peakAvg(jj));
   T2peakMin_poreSize(jj) = ((sqrt(Dtemp.*(2.*T2peakMinbParam.*rho^2+Dtemp)) - Dtemp)./rho)*10^3; 
   
   diffRegime(jj,1) = (rho*meanPoreSize(jj)*10^-3)./Dtemp;
   
   largeModeSum(jj,1) = sum(sliceT2(poreSizeDist > 2*10^-2));
   smallModeSum(jj,1) = sum(sliceT2(poreSizeDist < 2*10^-2));
   
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

   
   
   figure(4)
   hold on
   title('Exploring K estimate of pore size')
   grid on
   box on
   
   scatter(poreSizeDist,sliceT2,10,'filled')
   plot(meanPoreSize(jj),0.005,'r*','MarkerSize',10)
   
   plot(T2ML_poreSize(jj),0.005,'b*','MarkerSize',10)
   
   xlabel('Estimated pore size (mm)')
   ylabel('Amplitude')  
   set(gca,'FontSize',12)
   %plot(T2W_poreSize(jj),0.005,'g*','MarkerSize',10)
   %plot(T2peakAvg_poreSize(jj),0.005,'k*','MarkerSize',10)

   %plot(poreSizeDist(goodMinLoc),sliceT2(goodMinLoc),'b*','MarkerSize',10)

   
   %ind = find((widthRatio > mean(widthRatio,'omitNaN')) & (sortedPhi < mean(sortedPhi))')
  
   %plot(T2peakMax_poreSize(jj),0.005,'m*','MarkerSize',10)
   %plot(T2peakMin_poreSize(jj),0.005,'k*','MarkerSize',10)
   
   ylim([0,0.02])
    
   set(gca,'XScale','log')
   
   %pause
       
   close()
   
end
    
%%
%     for jj = 1:7
%         sliceT2 = sortedT2dist(jj,:);
% 
%         figure(5)
%         subplot(2,1,1)
%         box on 
%         grid on
%         hold on
%         %scatter(poreSizeDist, sliceT2,10,'r','filled')
%         plot(poreSizeDist,sliceT2,'r')
%         plot(meanPoreSize(jj),0.005,'r*','MarkerSize',10)
% 
%         ylim([0,0.02])
%         set(gca,'XScale','log')
%         
%     end
%     
%     for jj = 7:21
%         subplot(2,1,2)
%         sliceT2 = sortedT2dist(jj,:);
% 
%          box on 
%         grid on
%         hold on
%         %scatter(poreSizeDist, sliceT2,10,'b','filled')
%         plot(poreSizeDist,sliceT2,'b')
%         plot(meanPoreSize(jj),0.005,'r*','MarkerSize',10)
% 
%         
%         ylim([0,0.02])
%         set(gca,'XScale','log')
%         
%     end
%     
    figure(6)
    subplot(5,1,1)
    hold on
    title('Prominance Ratio')
    scatter(sortedK,promRatio,30,'filled')
    set(gca,'XScale','log')
    grid on
    box on
    
    subplot(5,1,2)
    hold on
    title('Width Ratio')
    scatter(sortedK,widthRatio,30,'filled')
    set(gca,'XScale','log')
    grid on
    box on
    
    subplot(5,1,3)
    hold on
    title('Porosity')
    scatter(sortedK,sortedPhi,30,'filled')
    set(gca,'XScale','log')
    grid on
    box on
    
    subplot(5,1,4)
    hold on
    
    title('Summed signal large mode')
    %scatter(sortedK, largeModeSum,30,'filled')
    %scatter(sortedK, smallModeSum, 30,'filled')
    
    scatter(sortedK, longSummedSignal,30,'filled')

    
    set(gca,'XScale','log')
    grid on
    box on
    
    subplot(5,1,5)
    hold on
    
    title('Summed signal small mode')
    scatter(sortedK, shortSummedSignal, 30,'filled')
    
    set(gca,'XScale','log')
    grid on
    box on
    %testData = promRatio((widthRatio > mean(widthRatio,'omitNaN')) & (sortedPhi < mean(sortedPhi))')
    %ind = find((widthRatio > mean(widthRatio,'omitNaN')) & (sortedPhi < mean(sortedPhi))')

    
    %%
    % Try to magically filter out the bimodal distributions that matter
    % Doesn't work very well, too many bimodal distributions that are
    % sensitive to the larger mode
%     
%     for jj = 1:length(T2estimator)
% 
%         sliceT2 = sortedT2dist(jj,:);
%        
%         
%         figure(7)
%         hold on
%         grid on
%         box on
%         
%         scatter(poreSizeDist,sliceT2,10,'filled')
%         
%         
%         plot(meanPoreSize(jj),0.005,'r*','MarkerSize',10)
%     
%         ylim([0,0.02])
%         set(gca,'XScale','log')
%         
%         
%         if (widthRatio(jj) > mean(widthRatio,'omitNan'))...
%                 && (promRatio(jj) > mean(promRatio,'omitNaN'))...
%                 && (smallModeSum(jj) > mean(smallModeSum))...
%                 && (largeModeSum(jj) < mean(largeModeSum))
%             
%             plot(T2peakMin_poreSize(jj),0.005,'k*','MarkerSize',10)
%         else
%            plot(T2W_poreSize(jj),0.005,'g*','MarkerSize',10)
%         end
%             
%         pause
%        
%         close()
%     end
% %         
    %%
%     % Try to figure out a range of pore sizes based on bimodal processing
%     for jj = 1:length(T2estimator)
%         
%         sliceT2 = sortedT2dist(jj,:);
% 
%         figure(7)
%         hold on
%         grid on
%         box on
%         
%         scatter(poreSizeDist,sliceT2,10,'filled')
%         plot(meanPoreSize(jj),0.005,'r*','MarkerSize',10)
%     
%         ylim([0,0.02])
%         set(gca,'XScale','log')
%         
%         if promRatio(jj) > 0.1
%             plot(T2peakMin_poreSize(jj),0.005,'k*','MarkerSize',10)
%             plot(T2W_poreSize(jj),0.005,'g*','MarkerSize',10)
%         else
%             plot(T2W_poreSize(jj),0.005,'g*','MarkerSize',10)
%         end
%     
%        % pause
%         close()
%     end
        
end


