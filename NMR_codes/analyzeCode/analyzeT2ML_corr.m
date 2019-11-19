% Correct T2ml data 

sites = {'wisc_all'}
legendNames = {'Wisc'}

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
    
    %%
    % Sort data and correct large T2ML
    
    [sortedLogT2ML, sortInd] = sort(logT2ML);
    sortedLogK = logK(sortInd,:);

    T2MLlogCutoff = -0.78;

    smallLogK = sortedLogK(sortedLogT2ML < T2MLlogCutoff);
    largeLogK = sortedLogK(sortedLogT2ML > T2MLlogCutoff);
    smallLogT2ML = sortedLogT2ML(sortedLogT2ML < T2MLlogCutoff);
    largeLogT2ML = sortedLogT2ML(sortedLogT2ML > T2MLlogCutoff);
    
    X = [ones(length(smallLogT2ML),1) smallLogT2ML];
    fitParams = X\(smallLogK);
           
    T2logSlope(1) = fitParams(2)
    T2logyInt(1) = fitParams(1);
    
    fitLogT2MLRange = logspace(-2,-0.8,200);
    fitLogT2MLRange = log10(fitLogT2MLRange);
    
    bestKFitSmall = T2logyInt(1) + (fitLogT2MLRange.*T2logSlope(1));
    
    % Try to correct y Int
    SDRyInt = (T2logyInt(1)/T2logSlope(1))*2;
    SDRestimate = SDRyInt + (fitLogT2MLRange.*2);
    
    figure(6)
    hold on
    
    %Try Seevers approx for T2
    T2B = 1;
    
    SeeversT2 = (T2ML.^(-1) - T2B^(-1)).^(-1);
    SeeversCheck = (T2ML.*T2B)./(T2B - T2ML);
    
    logSeeversT2 = log10(SeeversT2);
    
    scatter(logT2ML,logK,30,'filled')
    scatter(logSeeversT2, logK, 30,'filled')

%     scatter(T2ML,logK,30,'filled')
%     scatter(SeeversT2, logK, 30,'filled')
%     
    corrFactor = 0.4;
    
%      scatter(smallLogT2ML, smallLogK,30,'filled') 
%      scatter(largeLogT2ML - corrFactor, largeLogK, 30, 'Filled')
%      scatter(largeLogT2ML, largeLogK, 30, 'Filled')

    
    %plot(largeLogT2MLRange, bestKFitLarge,'LineWidth',2)
    %plot(smallLogT2MLRange, bestKFitSmall,'LineWidth',2)
    %plot(smallLogT2MLRange, SDRestimate,'k','LineWidth',2)

    
    grid on
    box on
    
    xlabel('log10 T2ML (s)')
    ylabel('log10 K (m/s)')
    
%     figure(7)
%     hold on
%     scatter(T2ML, K)
%     grid on
%     box on
%     xlabel(' T2ML (s)')
%     ylabel(' K (m/s)')
%     
    %%
    % Try fit on corrected data
    corrLargeT2ML = largeLogT2ML - corrFactor;
    corrlogT2ML = [smallLogT2ML; corrLargeT2ML];
    corrlogK = [smallLogK; largeLogK];
    
    X = [ones(length(corrlogT2ML),1) corrlogT2ML];
    fitParams = X\(corrlogK);
           
    T2logSlope(1) = fitParams(2)
    T2logyInt(1) = fitParams(1);
    
    fitLogT2MLRange = logspace(-2,-0.8,200);
    fitLogT2MLRange = log10(fitLogT2MLRange);
    
    bestKFit = T2logyInt(1) + (fitLogT2MLRange.*T2logSlope(1));
    
    figure(2)
    hold on
    
    scatter(corrlogT2ML, corrlogK,30,'Filled')
    plot(fitLogT2MLRange, bestKFit, 'k','LineWidth',2)
    
    grid on
    box on
    
    xlabel('log10 T2ML (s)')
    ylabel('log10 K (m/s)')
    
    xlim([-2, -0.2])
    
    %% Now calculate b
    corrK = 10.^(corrlogK);
    corrT2ML = 10.^(corrlogT2ML);
    
    corrb = SDR_b(corrK,m,n,phi,corrT2ML);
    
    figure(3)
    hold on
    
    scatter(corrK,corrb,30,'Filled')
    scatter(K,bProfile{kk},30,'Filled')
    
    grid on
    box on
    
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    
    xlabel('K (m/s)')
    ylabel('Corrected b')
    
    %% Plot T2ML vs K for both
    figure(4)
    subplot(121)
    
    hold on
    
    scatter(logT2ML, logK, 30,'Filled')
    
    grid on
    box on
    
    %set(gca,'XScale','log')
    %set(gca,'YScale','log')
    %set(gca,'XScale','log')
    
    ylabel('K (m/s)')
    xlabel('T2ML (s)')
    xlim([-1.5 -0.4])
    
    
    subplot(122)
    hold on
     
    X = [ones(length(corrlogT2ML),1) corrlogT2ML];
    fitParams = X\(corrlogK);
           
    T2logSlope(1) = fitParams(2)
    T2logyInt(1) = fitParams(1);
    
    fitLogT2MLRange = logspace(-2,-0.8,200);
    fitLogT2MLRange = log10(fitLogT2MLRange);
    
    bestKFit = T2logyInt(1) + (fitLogT2MLRange.*T2logSlope(1));
    
    scatter(corrlogT2ML, corrlogK, 30,'Filled')
    
    grid on
    box on
    
    %set(gca,'XScale','log')
    %set(gca,'YScale','log')
    
    ylabel('K (m/s)')
    xlabel('T2ML (s)')
    xlim([-1.5 -0.4])

    %%
    % Attempt to calculate mean pore radius from Godefroy et al 2001
    
    D = @(temp) (1.0413+0.039828*temp+0.00040318*temp^2)*10^-9;
    T2B = 2.2;
    alpha = 3;
    %bParam = (T2ML.*T2B.*alpha)./(T2B - T2ML);
    bParam = (corrT2ML.*T2B.*alpha)./(T2B - corrT2ML);
    rho = 44*10^-6; %m/s
    Dtemp = D(15);
    
    poreSize = (sqrt(Dtemp.*(2.*bParam.*rho^2+Dtemp)) - Dtemp)./rho;
    
    poreSizeLine = logspace(-7,-4,100);
    Kpred = (poreSizeLine.^4.5).*10^(19);
    
    figure(20)
    hold on
    scatter(poreSize, corrK,30,'Filled')
    plot(poreSizeLine, Kpred,'LineWidth',2)
    %scatter(poreSize, K,30,'Filled')
    
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    
    ylabel('K (m/s)')
    xlabel('Pore Size? (m)')
    
    ylim([10^-8,10^-2])
    xlim([10^-6,5*10^-5])
    
    grid on
    box on

    figure(21)
    hold on
    scatter(poreSize, corrT2ML,30,'Filled')

    set(gca,'XScale','log')
    set(gca,'YScale','log')

    ylabel('T2ML (s)')
    xlabel('Pore Size? (m)')

    grid on
    box on
    
end