% Analyze SDR b information
%close all

% sites = {'dpnmr_larned_east','dpnmr_larned_lwph','dpnmr_larned_west',...
%     'dpnmr_leque_east','dpnmr_leque_west','dpnmrA11','dpnmrA12',...
%     'dpnmrC1S','dpnmrC1SE','dpnmrC1SW'};
% 
% legendNames = {'Larned East','Larned lwph','Larned West','Leque East',...
%     'Leque West','A11','A12','C1S','C1SE','C1SW'};

% sites = {'adams_all','plainfield_all','dpnmr_larned_all','dpnmr_leque_all','gems_all'};
% legendNames = {'Adams','Plainfield','Larned','Leque','GEMS'};

%sites = {'Site1-WellG5','Site1-WellG6','Site2-WellPN1','Site2-WellPN2'}


sites = {'wisc_all','maurer_all'}
legendNames = {'Wisc','Maurer and Knight'}

%  sites = {'maurer_all'}
%  legendNames = {'Wisc'}

% All sites
% sites = {'dpnmr_larned_east','dpnmr_larned_lwph','dpnmr_larned_west',...
%     'dpnmr_leque_east','dpnmr_leque_west','dpnmrA11','dpnmrA12',...
%     'dpnmrC1S','dpnmrC1SE','dpnmrC1SW','Site1-WellG5','Site1-WellG6',...
%     'Site2-WellPN1','Site2-WellPN2'};

for kk = 1:length(sites)
    %close all
    baseName = sites{kk}
    baseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/USGS Data/';
    baseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/NMR-K-prediction/Data/Aggregated_Data/';

    %baseDir = 'I:\My Drive\USGS Project\USGS Data\';
    
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

    SDR_K = @(b,m,n,phi,T2ML) (b.*(phi.^m).*(T2ML).^n);
    SDR_b = @(K,m,n,phi,T2ML) K./((phi.^m).*(T2ML).^n);
    
    SOE_K = @(C,SOE,n) C.*(SOE.^n);
    SOE_C = @(K,SOE,n) K./(SOE.^n);
    
    m = 0;
    n = 2;
    
    bProfile{kk} = SDR_b(K,m,n,phi,T2ML);
    cProfile{kk} = SOE_C(K, SumEch, n);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Look at Hazen approximation
    hK = @(A,d) A.*d.^2; % in cm/s
    
    d10 = logspace(-3,1,1000); %in mm
    %d10 = 0.00001:0.0001:0.01 %in meters 0.1mm-1cm
    
    hazenK = hK(70,d10)*10^(-2); %in m/s   
    
    figure(23)

   % scatter(d10, hazenK, 40, 'filled')
    
    set(gca,'XScale','log')
    set(gca,'YScale','log')
        grid on
    box on
    
    xlabel('d10 (mm)')
    ylabel('K (m/s)')
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot b values versus hydraulic conducitivity 
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
    
    slope{kk} = fitParams(2)
    yInt{kk} = fitParams(1)
    
    KlinTest = logspace(-7, -2, 1000);
    KlogTest = log10(KlinTest);
    
    %KlinTest = 10.^(KlogTest);
    
    yInt{kk} = 0.4522;
    slope{kk} = 0.5522;
    
    bFitLineLog = (yInt{kk}) + KlogTest.*slope{kk};
    
    bFitLineLin = 10.^bFitLineLog;
    
    scatter(K, bProfile{kk},30,'Filled')
    %plot(KlinTest,bFitLineLin,'LineWidth',2,'HandleVisibility','off')
    
    
    %scatter(hazenK,d10, 40, 'filled')

    
    %plot(Krange, [minMaurer minMaurer],'b','LineWidth',2,'HandleVisibility','off')
    %plot(Krange, [maxMaurer maxMaurer],'b','LineWidth',2,'HandleVisibility','off')

    xlabel('Hydraulic Conductivity (m/s)')
    ylabel('Ideal SDR b Parameter (m/s^3)')
   

    set(gca,'xscale','log')
    set(gca,'yscale','log')
    set(gca,'FontSize',14)
    
%     xlim([10^-7 10^-2])
%     ylim([10^-5 10^-1])
%     
    

    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(21)
    subplot(121)
    grid on
    box on
    hold on
    
    scatter(z, bProfile{kk}, 40,'Filled')
    xlabel('Depth (m)')
    ylabel('Ideal SDR b Parameter')
    set(gca,'yscale','log')
    %ylim([10^-5 10^-1])
    set(gca,'FontSize',12)

    
    subplot(122)
    grid on
    box on
    hold on
    
    scatter(z, K, 40,'Filled')
    xlabel('Depth (m)')
    ylabel('Hydraulic Conductivity (m/s)')
    set(gca,'yscale','log')
    %ylim([10^-7 10^-2])
    set(gca,'FontSize',12)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Try filtering based on depth
    depthCutoff = 5;
    
    bValues = bProfile{kk};
    
    shallowK = K(z < depthCutoff);
    shallowb = bValues(z < depthCutoff);
    
    deepK = K(z > depthCutoff);
    deepb = bValues(z > depthCutoff);
    
    figure(22)
    hold on
    grid on
    box on
    
    %scatter(shallowK, shallowb,30,'Filled')
    %scatter(deepK, deepb,30,'Filled')
    scatter(cProfile{kk}, z, 30,'Filled')

    %plot(KlinTest,bFitLineLin,'LineWidth',2,'HandleVisibility','off')
    
    
    xlabel('Hydraulic Conductivity (m/s)')
    ylabel('Ideal SDR b Parameter')
   

    set(gca,'xscale','log')
    set(gca,'yscale','log')
    set(gca,'FontSize',16)
    
    %xlim([10^-7 10^-2])
    %ylim([10^-5 10^-1])
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     figure(3)
%     grid on
%     box on
%     hold on
%     
%     scatter(K, bProfile{kk},40,'Filled')
%     plot(KlinTest,bFitLineLin)
        
    
    figure(2)
%     minLine = plot(Krange, [minMaurer minMaurer],'b','LineWidth',2,'HandleVisibility','off');
%     maxLine = plot(Krange, [maxMaurer maxMaurer],'b','LineWidth',2,'HandleVisibility','off');
%    uistack(minLine,'bottom')
%    uistack(maxLine,'bottom')

    %legend([{'Maurer and Knight 2016 Range'}, {legendNames{kk}}],'Location','southeast')
    %legend(legendNames)

    % Try subtracting known b relation solved for previously
    bFitInterp = interp1(KlinTest,bFitLineLin,K);
    corrbProfile{kk} = bProfile{kk} - bFitInterp;
    
    figure(3)
    hold on
    grid on
    box on 
    
    scatter(K, corrbProfile{kk}, 30, 'Filled')
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

    h1 = histogram(log10(K),20)
    hold on
 
    xlabel('log K (m/s)')
    ylabel('Freq')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    % Let's try fitting K in two different regimes
    [sortedK, sortInd] = sort(K);
    sortedT2ML = T2ML(sortInd,:);
    
    %Kcutoff = 9*10^-5;
    Kcutoff = 10^-4;
    
    K_large = sortedK(sortedK > Kcutoff);
    K_small = sortedK(sortedK < Kcutoff);
    
    logK_small = log10(K_small);
    logK_large = log10(K_large);
    
    T2ML_large = sortedT2ML(sortedK > Kcutoff);
    T2ML_small = sortedT2ML(sortedK < Kcutoff);
    
    smallKrange = [10^-6 Kcutoff];
    largeKrange = [Kcutoff 10^-3];
    
    logSmallKrange = log10(smallKrange);
    logLargeKrange = log10(largeKrange);
    smallKFitLin = logspace(logSmallKrange(1),logSmallKrange(2),100);
    largeKFitLin = logspace(logLargeKrange(1),logLargeKrange(2),100);
    smallKFit = log10(smallKFitLin);
    largeKFit = log10(largeKFitLin);
    
    % Compute Linear Fit
    X = [ones(length(logK_small),1) logK_small];
    fitParams = X\T2ML_small;
           
    T2slope(1) = fitParams(2);
    T2yInt(1) = fitParams(1);
    
    X = [ones(length(logK_large),1) logK_large];
    fitParams = X\T2ML_large;
           
    T2slope(2) = fitParams(2);
    T2yInt(2) = fitParams(1);
    
   
    bestFitSmall = T2yInt(1) + (smallKFit.*T2slope(1));
    bestFitLarge = T2yInt(2) + (largeKFit.*T2slope(2));

    figure(28)
    subplot(121)
    hold on
    
    scatter(K_large, T2ML_large,40,'filled')
    scatter(K_small, T2ML_small,40,'filled')
    %scatter(sortedK, sortedT2ML,20,'filled')
    
    plot(smallKFitLin,bestFitSmall,'LineWidth',2)
    plot(largeKFitLin,bestFitLarge,'LineWidth',2)
    
    grid on
    box on
    set(gca,'xscale','log')
    ylabel('T2ML (s)')
    xlabel('K (m/s)')
    

    subplot(122)
    hold on
    
    
    X = [ones(length(logK_small),1) logK_small];
    fitParams = X\(T2ML_small.^2);
           
    T2slope(3) = fitParams(2);
    T2yInt(3) = fitParams(1);
    
    X = [ones(length(logK_large),1) logK_large];
    fitParams = X\(T2ML_large.^2);
           
    T2slope(4) = fitParams(2);
    T2yInt(4) = fitParams(1);
    
    bestFitSmall = T2yInt(3) + (smallKFit.*T2slope(3));
    bestFitLarge = T2yInt(4) + (largeKFit.*T2slope(4));
    
    scatter(K_large, T2ML_large.^2,40,'filled')
    scatter(K_small, T2ML_small.^2,40,'filled')
    %scatter(sortedK, sortedT2ML,20,'filled')
    
    plot(smallKFitLin,bestFitSmall,'LineWidth',2)
    plot(largeKFitLin,bestFitLarge,'LineWidth',2)
    
    grid on
    box on
    set(gca,'xscale','log')
    ylabel('T2ML (s)')
    xlabel('K (m/s)')
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    % Try choosing a "b" and then calculating K
    %
    figure()
    
    median_b = median(bProfile{kk});
    max_b = max(bProfile{kk});
    min_b = min(bProfile{kk});
    
    K_medianEst = SDR_K(median_b,m,n,phi,T2ML);
    K_maxEst = SDR_K(max_b,m,n,phi,T2ML);
    K_minEst = SDR_K(min_b,m,n,phi,T2ML);
    
    Kest = [K_medianEst K_maxEst K_minEst];
    Knames = [{'Median b'} {'Max b'} {'Min b'}];
    
    colors = {[0.9290,0.6940,0.1250],[0.466,0.6740,0.1880],[0.3,0.8,0.1]};

    plotKestKdpp_wdepth(K,Kest,Kest,z,Knames,colors)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Let's try to match K exactly
    figure()
    K_norm = SDR_K(0.0043,0,2,phi,T2ML);
    K_best = (K.^slope{kk}.*10.^yInt{kk}).*T2ML.^2;
    
    bestC = median(cProfile{kk});
    
    K_SOE = SOE_K(bestC,SumEch,2);
    
    K_est = [K_norm K_best];
    K_names = [{'K Normal b'},{'K Corrected b'},{'K SOE'}];
    plotKestKdpp_wdepth(K,K_est,K_est,z,K_names,...
        [0.5 0.5 0.5])
    
    K_factor = estimateKdiffFactor(K,K_est,1);
    K_factorAvg = mean(K_factor)
    
    legend(K_names)
    
    %%

% Filter the data by K    
%     [sortedLogK, sortInd] = sort(logK);
%     sortedLogT2ML = logT2ML(sortInd,:);
%     
%     %KlogCutoff = -4.3;    
%     KlogCutoff = -3
%     
%     largeLogK = sortedLogK(sortedLogK > KlogCutoff);
%     largeLogT2ML = sortedLogT2ML(sortedLogK > KlogCutoff);
%     
%     smallLogK = sortedLogK(sortedLogK < KlogCutoff);
%     smallLogT2ML = sortedLogT2ML(sortedLogK < KlogCutoff);
%     
%     X = [ones(length(largeLogT2ML),1) largeLogT2ML];
%     fitParams = X\(largeLogK);
%            
%     T2logSlope(1) = fitParams(2);
%     T2logyInt(1) = fitParams(1);
%     
%     X = [ones(length(smallLogT2ML),1) smallLogT2ML];
%     fitParams = X\(smallLogK);
%            
%     T2logSlope(2) = fitParams(2);
%     T2logyInt(2) = fitParams(1);
%     
%     
%     largeLogT2MLRange = logspace(-1.4,-0.4,100);
%     largeLogT2MLRange = log10(largeLogT2MLRange);
%     
%     smallLogT2MLRange = logspace(-2,-1,100);
%     smallLogT2MLRange = log10(smallLogT2MLRange);
%     
%     bestKFitLarge = T2logyInt(1) + (largeLogT2MLRange.*T2logSlope(1));
%     bestKFitSmall = T2logyInt(2) + (smallLogT2MLRange.*T2logSlope(2));

% Filter the data by T2ML
    %close all
    
    [sortedLogT2ML, sortInd] = sort(logT2ML);
    sortedLogK = logK(sortInd,:);

    T2MLlogCutoff = -0.9;

    smallLogK = sortedLogK(sortedLogT2ML < T2MLlogCutoff);
    largeLogK = sortedLogK(sortedLogT2ML > T2MLlogCutoff);
    smallLogT2ML = sortedLogT2ML(sortedLogT2ML < T2MLlogCutoff);
    largeLogT2ML = sortedLogT2ML(sortedLogT2ML > T2MLlogCutoff);
    
    X = [ones(length(smallLogT2ML),1) smallLogT2ML];
    fitParams = X\(smallLogK);
           
    T2logSlope(1) = fitParams(2)
    T2logyInt(1) = fitParams(1);
    
    smallLogT2MLRange = logspace(-2,-0.8,200);
    smallLogT2MLRange = log10(smallLogT2MLRange);
    
    bestKFitSmall = T2logyInt(1) + (smallLogT2MLRange.*T2logSlope(1));
    
    % Try to correct y Int
    SDRyInt = (T2logyInt(1)/T2logSlope(1))*2;
    SDRestimate = SDRyInt + (smallLogT2MLRange.*2);
    
    figure(6)
    hold on
    
    %Try Seevers approx for T2
    T2B = 2;
    SeeversT2 = (T2ML.^(-1) - T2B^(-1)).^(-1);
    logSeeversT2 = log10(SeeversT2);
    
    %scatter(logT2ML,logK,30,'filled')
    %scatter(logSeeversT2, logK, 30,'filled')
    
    corrFactor = 0.5;
    
    scatter(smallLogT2ML, smallLogK,30,'filled') 
    scatter(largeLogT2ML - corrFactor, largeLogK, 30, 'Filled')
    %scatter(largeLogT2ML, largeLogK, 30, 'Filled')

    
    %plot(largeLogT2MLRange, bestKFitLarge,'LineWidth',2)
    plot(smallLogT2MLRange, bestKFitSmall,'LineWidth',2)
    plot(smallLogT2MLRange, SDRestimate,'k','LineWidth',2)
    

    
    grid on
    box on
    
    xlabel('log10 T2ML (s)')
    ylabel('log10 K (m/s)')
    
    
    
end

% 
% figure(2)
% minLine = plot(Krange, [minMaurer minMaurer],'b','LineWidth',2,'HandleVisibility','off');
% maxLine = plot(Krange, [maxMaurer maxMaurer],'b','LineWidth',2);
% uistack(minLine,'bottom')
% uistack(maxLine,'bottom')
% 
% legend(['Maurer and Knight 2016 Range', legendNames])