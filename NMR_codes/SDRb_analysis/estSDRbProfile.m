% Compare b estimates from Kansas/Washington data

% sites = {'dpnmr_larned_east','dpnmr_larned_lwph','dpnmr_larned_west',...
%     'dpnmr_leque_east','dpnmr_leque_west','dpnmrA11','dpnmrA12',...
%     'dpnmrC1S','dpnmrC1SE','dpnmrC1SW'};
% 
% legendNames = {'Larned East','Larned lwph','Larned West','Leque East',...
%     'Leque West','A11','A12','C1S','C1SE','C1SW'};

sites = {'adams_all','plainfield_all','dpnmr_larned_all','dpnmr_leque_all','gems_all'};
legendNames = {'Adams','Plainfield','Larned','Leque','GEMS'};
% 
% sites = {'wisc_all'}
% legendNames = {'wisc'}

for kk = 1:length(sites)
    close all
    baseName = sites{kk}
    %baseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/USGS Data/';
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

    nmrName
    [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(nmrName);
    
    SDR_b = @(K,m,n,phi,T2ML) K./((phi.^m).*(T2ML).^n);

    m = 0;
    n = 2;
    
    bProfile{kk} = SDR_b(K,m,n,phi,T2ML);
    
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
    
    slope{kk} = fitParams(2)
    yInt{kk} = fitParams(1)
    
    KlinTest = logspace(-6, -3, 100);
    KlogTest = log10(KlinTest);
    
    %KlinTest = 10.^(KlogTest);
    
    bFitLineLog = (yInt{kk}) + KlogTest.*slope{kk};
    
    bFitLineLin = 10.^bFitLineLog;
    
    scatter(K, bProfile{kk},40,'Filled')
    plot(KlinTest,bFitLineLin, 'HandleVisibility','off')
    
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
    maxLine = plot(Krange, [maxMaurer maxMaurer],'b','LineWidth',2);
    uistack(minLine,'bottom')
    uistack(maxLine,'bottom')

    legend([{'Maurer and Knight 2016 Range'}, {legendNames{kk}}],'Location','southeast')

    fileName = strcat(baseName,'_bRange.fig');
    savefig(fileName)
    
    fileName = strcat(baseName,'_bRange.png');
    print('-dpng','-r300',fileName)
    
end

% 
% figure(2)
% minLine = plot(Krange, [minMaurer minMaurer],'b','LineWidth',2,'HandleVisibility','off');
% maxLine = plot(Krange, [maxMaurer maxMaurer],'b','LineWidth',2);
% uistack(minLine,'bottom')
% uistack(maxLine,'bottom')
% 
% legend(['Maurer and Knight 2016 Range', legendNames])