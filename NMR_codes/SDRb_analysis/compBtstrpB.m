% Bootstrap b comp to ideal SDR b

%% Compare to bstrpB for m = 0, n = 2
% sites = {'Site1-WellG5','Site1-WellG6','Site2-WellPN1','Site2-WellPN2'}
% bstrpB = [0.0040 0.0024 0.0075 0.0047]
% color = [0 114 189]/255;
% colorStar = [198 158 190]/255;
% 

sites = {'dpnmr_larned_east','dpnmr_larned_lwph','dpnmr_larned_west',...
    'dpnmr_leque_east','dpnmr_leque_west','dpnmrA11','dpnmrA12',...
    'dpnmrC1S','dpnmrC1SE','dpnmrC1SW'};

bstrpB = [0.027,0.024,0.016,0.027,0.022,0.032,0.032,0.028,0.028,0.028];
color = [217 83 25]/255;
colorStar = [196 209 154]/255;

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
    n = 1;
    
    bProfile{kk} = SDR_b(K,m,n,phi,T2ML);
    cProfile{kk} = SOE_C(K, SumEch, n);
    T2mlProfile{kk} = T2ML;
    
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
    
    logMeanK = 10.^mean(log10(K));
    logMeanT2ML = 10.^mean(log10(T2ML));
    
    scatter(K, bProfile{kk},25,'MarkerEdgeColor',color, 'MarkerFaceColor',color)
    plot(logMeanK,bstrpB(kk),'s','MarkerSize',15,'MarkerEdgeColor','k','LineWidth',2)
    
    %plot(logMeanK,bstrpB(kk),'p','MarkerSize',20,'MarkerFaceColor',colorStar,'MarkerEdgeColor','k')

    %plot(KlinTest,bFitLineLin,'LineWidth',2,'HandleVisibility','off')
    
    
    %scatter(hazenK,d10, 40, 'filled')

    
    %plot(Krange, [minMaurer minMaurer],'b','LineWidth',2,'HandleVisibility','off')
    %plot(Krange, [maxMaurer maxMaurer],'b','LineWidth',2,'HandleVisibility','off')

    xlabel('Hydraulic Conductivity (m/s)')
    ylabel('Ideal SDR b_K Parameter (m/s^3)')
   

    set(gca,'xscale','log')
    set(gca,'yscale','log')
    set(gca,'FontSize',14)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot b vs T2ML 
    
    figure(3)
    box on 
    grid on
    hold on
    
    minMaurer = 1.5E-02;
    maxMaurer = 3.6E-02;
    Krange = [10^-6 10^-3];
    
    scatter(T2mlProfile{kk}, bProfile{kk},25,'MarkerEdgeColor',color, 'MarkerFaceColor',color)
    plot(logMeanT2ML,bstrpB(kk),'p','MarkerSize',20,'MarkerFaceColor',colorStar,'MarkerEdgeColor','k')
    
    xlabel('T2ML (s)')
    ylabel('Ideal SDR b_K Parameter (m/s^3)')
    
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    set(gca,'FontSize',14)
%     xlim([10^-7 10^-2])
%     ylim([10^-5 10^-1])

end
%     