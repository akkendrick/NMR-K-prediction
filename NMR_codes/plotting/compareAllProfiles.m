% Compare profiles to one another 

clear
close all

load bProfileData.mat
load bestKModels_old.mat

sites = {'Site1-WellG5','Site1-WellG6',...
    'Site2-WellPN1','Site2-WellPN2'};
 
%waterTable = [2.0469,2.1248,5.0285,4.7476]; % rel ground surface
%waterTable = [3.004,2.963,5.727,5.408]; % rel top of casing
depthOffsets = [0.75,0.95,0.75,0.75];
waterTable = [1.181,1.181,5.0285,4.7476];

for kk = 1:length(sites)
    baseName = sites{kk}
    baseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/USGS Data/';
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
        site = 'Site1-WellG6';
        name = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';
        nmrName = name;
    elseif strcmp(baseName,'Site1-WellG6above')
        site = 'Site1-WellG6';
        name = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';
        nmrName = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1_above'
    elseif strcmp(baseName,'Site1-WellG6below')
        site = 'Site1-WellG6';
        name = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';
        nmrName = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1_below'
    elseif strcmp(baseName,'Site2-WellPN1')
        site = 'Site2-WellPN1';
        name = 'Pl_W1_Tr5_20x_MPp75aLS_F1n2_wRIN_wRFI_Reg50_Va1';
        nmrName = name;
    elseif strcmp(baseName,'Site2-WellPN2')
        site = 'Site2-WellPN2';
        name = 'W2_Tr5_20x_MPp75aLS_Reg50_wRIN_wRFI_Va1';
        nmrName = name;
    end

    in1 = [baseDir site '/' name '/' name '_SE_decay' '.txt']; 
    in2 = [baseDir site '/' name '/' name '_1Dvectors.txt'];
    in3 = [baseDir site '/' name '/' strcat(site,'_DPP_filt.txt')];
    in4 = [baseDir site '/' name '/' name '_SE_decay_time.txt'];
    in5 = [baseDir site '/' name '/' name '_T2_dist.txt'];
    in6 = [baseDir site '/' name '/' name '_T2_bins_log10s.txt'];

    load('sonicCoreT2B.mat','sonicCoreT2BData')
    T2B_depth = sonicCoreT2BData.Depthm;
    T2B_depth = flipud(T2B_depth);  
    
    T2B_peak = sonicCoreT2BData.T2Bpeak*(10^-3);
    T2B_peak = flipud(T2B_peak); 
    
    decaycurve = load(in1); 
    dparam = dlmread(in2,'\t',1,0); 
    DPPdat = load(in3); 
    t = load(in4);
    
    T2dist = load(in5);
    T2dist(:,1) = T2dist(:,1) - depthOffsets(kk);
    
    T2logbins = load(in6);
    % set needed variables
    NMRphi = dparam(:,2);

    S = decaycurve(:,2:end); 

    %Dk = olddat(:,4)*1.16e-5; 
    %z_dk = olddat(:,1); 
    
    %colors 

    [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(nmrName); 

    Dk = DPPdat(:,2)*1.16e-5; % converts K from m/day to m/s
    z_dk = DPPdat(:,1);

    origT2dist = T2dist;
    depths = T2dist(:,1);
    
    waterTableLine = [waterTable(kk) waterTable(kk)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot NMR porosity
    figure(1)
    subplot(151)
    box on
    grid on

    hold on
    plot(NMRphi,depths,'LineWidth',2)
    %scatter(NMRphi,depths,40,'Filled')
    %plot([0,1],waterTableLine,'r','LineWidth',2)
     
    xlabel('NMR Porosity')
    set(gca, 'YDir','reverse')
    ylim([min(depths),max(depths)])
    xlim([0,0.6])
    set(gca,'FontSize',14)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot T2ML
    subplot(152)    
    waterTableLine = [waterTable(kk) waterTable(kk)];
    
    box on
    grid on
    hold on
    %plot([-1.5,0],waterTableLine,'r','LineWidth',3,'HandleVisibility','off')

    scatter(log10(T2ML),z,60,'Filled')
    
    set(gca,'FontSize',14)
    set(gca, 'YDir','reverse')
    ylim([min(depths),max(depths)])
    xlabel('log_{10} T_{2ML}')
    xlim([-1.5 0])
    set(gca,'YTickLabel',[])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot b profile
    subplot(153)
    grid on
    box on
    hold on
    waterTableLine = [waterTable(kk) waterTable(kk)];
    
    %plot([0,1],waterTableLine,'r','LineWidth',3,'HandleVisibility','off')
    scatter((bProfile{kk}),z,60,'Filled')
     
    set(gca,'YDir','reverse')
    %xlim([min(bProfile),max(bProfile)])
    ylim([min(depths),max(depths)])
    xlim([10^-4 10^0])
    set(gca,'XScale','log')
    xlabel('SDR b')
    set(gca,'FontSize',14)
    set(gca,'YTickLabel',[])

    disp('Mean/Min/Max b value')
    mean(bProfile{kk})

    bProfileTot{kk,:} = bProfile;
    
    %title(baseName)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot hydraulic conductivity
    subplot(154)
    box on
    grid on
    
    hold on
    scatter(K,z,60,'Filled')
    %scatter(kProfile{:,kk},z,60,'Filled')
    
    xlabel('Hydraulic Conductivity (m/s)')
    set(gca, 'YDir','reverse')
    ylim([min(depths),max(depths)])
    set(gca,'XScale','log')
    set(gca,'FontSize',14)
    xlim([10^-6 10^-3])
    set(gca,'YTickLabel',[])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot error profile
    subplot(155)
    box on
    grid on
    hold on
    
    
    error = abs(bestSDRModel{kk} - K)
    %kFactor = bestSDRModel{kk} ./ K;
    
    scatter(error,z,60,'Filled')
    %scatter(kFactor,z,60,'Filled')
    
    xlabel('Error')
    set(gca, 'YDir','reverse')
    ylim([min(depths),max(depths)])
    set(gca,'FontSize',14)
    %xlim([10^-5 10^-3])
    set(gca,'YTickLabel',[])
    set(gca,'XMinorTick','on')
    %set(gca,'XScale','log')

    
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
    
    KlinTest = logspace(-6, -3, 100);
    KlogTest = log10(KlinTest);
    
    %KlinTest = 10.^(KlogTest);
    
    bFitLineLog = (yInt{kk}) + KlogTest.*slope{kk};
    
    bFitLineLin = 10.^bFitLineLog;
    
    scatter(K, bProfile{kk},40,'Filled')
    %plot(KlinTest,bFitLineLin)
    
    %plot(Krange, [minMaurer minMaurer],'b','LineWidth',2,'HandleVisibility','off')
    %plot(Krange, [maxMaurer maxMaurer],'b','LineWidth',2,'HandleVisibility','off')

    xlabel('Hydraulic Conductivity (m/s)')
    ylabel('Ideal SDR b Parameter')
   
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    set(gca,'FontSize',14)
    
    figure(3)
    grid on
    box on
    hold on
    
    scatter(K, bProfile{kk},40,'Filled')
    plot(KlinTest,bFitLineLin)
    
    disp('test')
    

end

figure (1)
subplot(151)
grid minor

subplot(155)
legend(sites)
grid minor

figure(2)
minLine = plot(Krange, [minMaurer minMaurer],'b','LineWidth',2,'HandleVisibility','off');
maxLine = plot(Krange, [maxMaurer maxMaurer],'b','LineWidth',2);
uistack(minLine,'bottom')
uistack(maxLine,'bottom')

legend(['Maurer and Knight 2016 Range', sites])
