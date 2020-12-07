% Analyze SDR b information
clear
%close all


sites = {'Site1-WellG5','Site1-WellG6','Site2-WellPN1','Site2-WellPN2'};

%sites = {'Site1-WellG5'};
%legendNames = [{'Site1-WellG5'}];

waterTable = [2.1248,2.0469,5.0285,4.7476]; % rel ground surface
depthOffsets = [0.95,0.75,0.75,0.75];

for kk = 1:length(sites)
    %close all
    baseName = sites{kk}
    
    rawBaseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/USGS Data/';
    baseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/NMR-K-prediction/Data/Aggregated_Data/';
    
    
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
        nmrName = baseName;
    end

    [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(nmrName);

    in1 = [rawBaseDir site '/' name '/' name '_T2_dist' '.txt']; 
    in2 = [rawBaseDir site '/' name '/' name '_T2_bins_log10s' '.txt']; 
    in3 = [rawBaseDir site '/' name '/' name '_1Dvectors' '.txt'];

    T2dist = load(in1);
    T2logbins = load(in2);
    T2linbins = 10.^(T2logbins);


    T2dist(:,1) = T2dist(:,1) - depthOffsets(1); 
       
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
       %filtPhi(jj) = phi(index)
    end
    
    [sortedK, sortInd] = sort(K);
    sortedT2dist = filtT2dist(sortInd,:);
    sortedPhi = phi(sortInd);
    
    
    figure()
    ribbon(sortedT2dist')
    title(baseName)
    
    xTicks = [0 2 4 6 8 10 12 14 16 18 20];
    yTicks = [0 10 20 30 40 50 60 70 80 90 100];
    
    newXLabel = [0; sortedK(2:2:end)];
    newYLabel = [0 T2axis(10:10:end)];
    
    % Try plotting large and small K separately
    figure(31)
    subplot(2,1,1)
    highK = sortedK(sortedK > 10^-4);
    highT2dist = sortedT2dist(sortedK > 10^-4,:);
    
    ribbon(highT2dist')
        title('High K')

    subplot(2,1,2)
    lowK = sortedK(sortedK < 10^-4);
    lowT2dist = sortedT2dist(sortedK < 10^-4,:);
    
    ribbon(lowT2dist')
    title('Low K')

    
%     set(gca,'XTick',xTicks)
%     set(gca,'YTick',yTicks)
%     
%     set(gca,'XTickLabel',newXLabel)
%     set(gca,'YTickLabel',newYLabel)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Try pulling portion of T2 distribution smaller than 60 ms (60.0350)
    
    T2bin = 51; 
    shortT2 = sortedT2dist(:,1:T2bin);
    longT2 = sortedT2dist(:,T2bin:end);
    
    figure(50)
    ribbon(shortT2')
    
    smoothFunc = @(x,a,b) 0.5+0.5*tanh((x-a)/b);

    x = 1:1:T2bin;
    weighting = fliplr(smoothFunc(x,25,14))*40;
    
    weightedShortT2 = shortT2.*weighting;
    
    figure(55)
    ribbon(weightedShortT2')
    
    shortT2Factor = sum(weightedShortT2,2);
    
    
    for kk = 1:length(K)
        shortT2ML_profile(kk) = sum(T2data(filtIndices(kk),1:T2bin).*(T2logbins(1:T2bin)))./sum(T2data(filtIndices(kk),1:T2bin));
        %T2ML_profile(kk) = sum(T2data(filtIndices(kk),:).*((T2logbins)))./sum(T2data(filtIndices(kk),:));
    end
    
    for kk = 1:length(T2depths)
        T2ML_profile(kk) = sum(T2data(kk,:).*((T2logbins)))./sum(T2data(kk,:));
    end
    
    figure(51)
    subplot(3,1,1)
    
    scatter(flipud(K),10.^shortT2ML_profile,40,'Filled')
    xlabel('K (m/s)')
    ylabel('T_{2ML} of short comp')
    title(baseName)

    
    set(gca,'XScale','log')
    
    grid on
    box on
    
    subplot(3,1,2)
    hold on
    
    summedShortSignal = sum(shortT2,2);
    summedLongSignal = sum(longT2,2);
    
    scatter(flipud(K),summedShortSignal,40,'Filled');
    scatter(flipud(K),summedLongSignal,40,'Filled');
    
    ylabel('Summed signal amplitude')
    xlabel('K (m/s)')
    legend({'Short Signal','Long Signal'})
    
    grid on 
    box on
    set(gca,'XScale','log')
    set(gca,'YScale','log')

    subplot(3,1,3)
    hold on
    
    grid on
    box on 
    %scatter(K,T2ML,40,'Filled')
    scatter((K), 10.^T2ML_profile(filtIndices),40,'Filled')
    set(gca,'XScale','log')
    ylabel('T_{2ML} (s)')
    xlabel('K (m/s)')
    
    figure(52)
    for kk = 1:length(K)
    
        plot((T2logbins(1:T2bin)), shortT2(kk,:))
        
        hold on

        scatter(shortT2ML_profile(kk),4*10^-3)
        
        hold off
        %disp('Yah')
    end
    %%
    % Quickly estimate T2 min/max estimator

    % Attempt to estimate pore radius from Godefroy
    D = @(temp) (1.0413+0.039828*temp+0.00040318*temp^2)*10^-9;
    T2B = 2.2;
    alpha = 3;

    rho = 44*10^-6; %m/s
    rhoMod = 10 ;
    rho = rho * rhoMod;
    
    Dtemp = D(15);
    
    
    for jj = 1:length(z)
        sliceT2 = sortedT2dist(jj,:);
        T2W(jj,1) = sum(T2linbins.*sliceT2)./phi(jj);
    end    
    
    T2estimator = T2W;
    for jj = 1:length(T2estimator)
        sliceT2 = sortedT2dist(jj,:);
        
        bParam = (T2linbins.*T2B.*alpha)./(T2B - T2linbins);
        poreSizeDist = ((sqrt(Dtemp.*(2.*bParam.*rho^2+Dtemp)) - Dtemp)./rho)*10^3;    
        
        figure()
        plot(poreSizeDist)
        
    % Identify bimodal distributions
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
       
       largeModeSum(jj,1) = sum(sliceT2(poreSizeDist > 2*10^-2));
       smallModeSum(jj,1) = sum(sliceT2(poreSizeDist < 2*10^-2));
       
    end

     % Try to figure out a range of pore sizes based on bimodal processing
   
    for jj = 1:length(T2estimator)
        if (widthRatio(jj) > mean(widthRatio,'omitNan'))...
                && (promRatio(jj) > mean(promRatio,'omitNaN'))...
                && (smallModeSum(jj) > mean(smallModeSum))...
                && (largeModeSum(jj) < mean(largeModeSum))            
            
            T2_lowerLim(jj,1) = T2peakMin(jj);
            T2_upperLim(jj,1) = T2W(jj);
        else
            T2_lowerLim(jj,1) = T2W(jj);
            T2_upperLim(jj,1) = T2W(jj);
        end
    end

    
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Try plotting SumEch
    
    figure(53)

    
    scatter(K,SumEch,40,'Filled')
    xlabel('K (m/s)')
    ylabel('Sum of Echoes')
        grid on
    box on
    set(gca,'XScale','log')
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
       
    dparam = dlmread(in3,'\t',1,0);   
    
    SDR_b = @(K,m,n,phi,T2ML) K./((phi.^m).*(T2ML).^n);

    m = 0;
    n = 2;
    
    %bProfile{kk} = SDR_b(K,m,n,phi,T2ML);
    bProfile{kk} = SDR_b(sortedK,m,n,sortedPhi,T2_lowerLim);

    
    
    % Plot b values versus hyraulic conducitivity 
    figure(2)
    box on 
    grid on
    hold on
    
    minMaurer = 1.5E-02;
    maxMaurer = 3.6E-02;
    Krange = [10^-6 10^-3];
    
    KLog = log10(K);
    bLog = log10(bProfile{kk});
    
    % Compute linear fit
    X = [ones(length(KLog),1) KLog];
    fitParams = X\bLog;
    
    slope{kk} = fitParams(2);
    yInt{kk} = fitParams(1);
    
    KlinTest = logspace(-7, -2, 1000);
    KlogTest = log10(KlinTest);
    
    %KlinTest = 10.^(KlogTest);
    
    bFitLineLog = (yInt{kk}) + KlogTest.*slope{kk};
    
    bFitLineLin = 10.^bFitLineLog;
    
    scatter(K, bProfile{kk}',40,'Filled')
    plot(KlinTest,bFitLineLin,'LineWidth',2,'HandleVisibility','off')
    
    %plot(Krange, [minMaurer minMaurer],'b','LineWidth',2,'HandleVisibility','off')
    %plot(Krange, [maxMaurer maxMaurer],'b','LineWidth',2,'HandleVisibility','off')

    xlabel('Hydraulic Conductivity (m/s)')
    ylabel('Ideal SDR b Parameter')
   
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    set(gca,'FontSize',14)
    
    xlim([10^-7 10^-2])
    ylim([10^-5 10^1])
    
%     figure(3)
%     grid on
%     box on
%     hold on
%     
%     scatter(K, bProfile{kk},40,'Filled')
%     plot(KlinTest,bFitLineLin)      
    
    figure(2)
    minLine = plot(Krange, [minMaurer minMaurer],'b','LineWidth',2,'HandleVisibility','off');
    maxLine = plot(Krange, [maxMaurer maxMaurer],'b','LineWidth',2,'HandleVisibility','off');
    uistack(minLine,'bottom')
    uistack(maxLine,'bottom')

    %legend([{'Maurer and Knight 2016 Range'}, {legendNames{kk}}],'Location','southeast')
    legend(legendNames)

    % Try subtracting known b relation solved for previously
    bFitInterp = interp1(KlinTest,bFitLineLin,K);
    corrbProfile{kk} = bProfile{kk} - bFitInterp;
    
    figure(3)
    hold on
    grid on
    box on 
    
    scatter(K, corrbProfile{kk}, 40, 'Filled')
    meanCorrb(kk) = mean(corrbProfile{kk});
    stdCorrb(kk) = std(corrbProfile{kk});
    
    figure(4)
    subplot(231)
    hold on
    box on
    grid on
    
    scatter(K,phi,40,'filled')
    set(gca,'xscale','log')
    xlabel('K (m/s)')
    ylabel('\phi')
    
    subplot(232)
    hold on 
    grid on
    box on
    
    scatter(K,z,40,'filled')
    set(gca,'xscale','log')
    ylabel('z (m)')
    xlabel('K (m/s)')
    
    subplot(233)
    hold on 
    grid on
    box on
    
    scatter(K,SumEch,40,'filled')
    set(gca,'xscale','log')
    ylabel('Sum of Echoes')
    xlabel('K (m/s)')   
    
    subplot(234)
    hold on 
    grid on
    box on
    
    scatter(K,T2ML,40,'filled')
    set(gca,'xscale','log')
    ylabel('T2ML (ms)')
    xlabel('K (m/s)')   
    
    subplot(235)
    hold on
    grid on
    box on
    
    scatter(K,T2ML.^2,40,'filled')
    set(gca,'xscale','log')
    ylabel('T2ML^2 (ms^2)')
    xlabel('K (m/s)') 
    
    subplot(236)
    hold on
    grid on
    box on
    scatter(sortedK, shortT2Factor,40','filled')
    set(gca,'xscale','log')
    ylabel('Short T_2 Factor')
    xlabel('K (m/s)') 
    
    
end

% 
% figure(2)
% minLine = plot(Krange, [minMaurer minMaurer],'b','LineWidth',2,'HandleVisibility','off');
% maxLine = plot(Krange, [maxMaurer maxMaurer],'b','LineWidth',2);
% uistack(minLine,'bottom')
% uistack(maxLine,'bottom')
% 
% legend(['Maurer and Knight 2016 Range', legendNames])