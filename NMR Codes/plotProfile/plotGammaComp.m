 sites = {'Site2-WellPN2','Site2-WellPN2'}
 
gammaNames = {'Site1-G5-Well1-gamma-EMI-bLS.csv','Site1-G6-Well2-gamma-EMI-bLS.csv'};
    
gammaNames = {'Site1-G6-Well2-gamma-EMI-bLS.csv','Site2-well2_EMI and Gamma_MbTOC.csv'};

waterTable = [2.0469,2.1248,5.0285,4.7476]; % rel ground surface
%waterTable = [3.004,2.963,5.727,5.408]; % rel top of casing
depthOffsets = [0.75,0.95,0.75,0.75];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 for kk = 1:length(sites)
    baseName = sites{kk}
%      baseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/USGS Data/';
%      gammaBaseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/USGS Data/WI_gamma-EMI-bLS_csvfiles';

    baseDir = 'I:\My Drive\Stanford\USGS Project\Field Data\USGS Data\';
    gammaBaseDir = 'I:\My Drive\Stanford\USGS Project\Field Data\USGS Data\WI_gamma-EMI-bLS_csvfiles';

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
    in3 = [baseDir site '/' name '/' strcat(site,'_DPP.txt')];
    in4 = [baseDir site '/' name '/' name '_SE_decay_time.txt'];
    in5 = [baseDir site '/' name '/' name '_T2_dist.txt'];
    in6 = [baseDir site '/' name '/' name '_T2_bins_log10s.txt'];
    in7 = [gammaBaseDir '/' gammaNames{kk}];
    
    gammaEMIdata = csvread(in7,3,0);

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

    [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(nmrName); 

    Dk = DPPdat(:,2)*1.16e-5; % converts K from m/day to m/s
    z_dk = DPPdat(:,1);

    origT2dist = T2dist;
    
    %%
    zNMR = z';

    figure(1)
    hold on

    grid on 
    box on
    
    gammaEMIdepth = gammaEMIdata(:,1);
    gamma = gammaEMIdata(:,2);
    EMI = gammaEMIdata(:,3);

    EMIdepth = gammaEMIdepth(EMI > 0);
    goodEMI = EMI(EMI > 0);

    %     origPosition = get(0,'defaultfigureposition');
    %     ax1 = axes('Position',origPosition);

    smoothWindow = 10;
    smoothGamma = smooth(gamma, smoothWindow);

    %plot(gamma, depth, 'parent', ax1)
    plot(smoothGamma, gammaEMIdepth,'LineWidth',2)

    %ylim([min(z),max(z)])
    xlim([0,100])

    set(gca,'YDir','reverse')


    xlabel('Counts/second')
    ylabel('Depth (m)')
    set(gca,'FontSize',12)
    
    figure(2)
    hold on 
    plot(goodEMI, EMIdepth,'LineWidth',2)
    set(gca,'XScale','log')
    set(gca,'YDir','reverse')
    grid on 
    box on
    
    
 end